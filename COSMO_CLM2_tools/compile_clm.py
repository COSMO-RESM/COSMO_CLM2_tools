from argparse import ArgumentParser, RawTextHelpFormatter
from glob import glob
from subprocess import check_call
import os
from shutil import rmtree


def compile_clm():

    # Define and parse command line arguments
    # ---------------------------------------

    dsc = "Compile CLM on Piz Daint. A case will be created in a subfolder of your ${SCRATCH}.\n"\
          " WARNING: tool has to be run from the default Prg-Env-cray environment"
    parser = ArgumentParser(description=dsc, formatter_class=RawTextHelpFormatter)
    parser.add_argument('cesm_trunk', help="path to the CESM directory")
    parser.add_argument('--clm_version', choices=['4.0', '4.5'], default='4.0', help="CLM version")
    parser.add_argument('-c', '--compiler', help="compiler to use (default: pgi)", default='pgi')
    parser.add_argument('-v', '--compiler_version', help="switch to this version of the compiler\n"\
                        "This is not recommended by CSCS")
    parser.add_argument('-d', '--debug', help="compile in debug mode (default: false)",
                        action='store_true')
    parser.add_argument('--src_mod', action='append',
                        help="path to additionnal/modified sources (e.g. oasis interface)\n"\
                        "has to be a folder containing src.xxx subfolders, can be specified several times")
    parser.add_argument('-o', '--output', help="output executable file path (default: ./cesm.exe)",
                        default='./cesm.exe')
    parser.add_argument('--no_exe', help="do not execute build_cesm.bash, leave it to any suited modification before actual compilation.",
                        action='store_false', dest='execute')
    opts = parser.parse_args()


    # Init some variables
    # -------------------

    CESM_TRUNK = opts.cesm_trunk
    EXP = 'clm{:s}_bld'.format(opts.clm_version)
    CASEDIR = os.path.join(os.environ['SCRATCH'], EXP)
    if os.path.exists(CASEDIR):
        rmtree(CASEDIR)
    RES = '1.9x2.5_gx1v6'
    COMP = 'ITEST'
    MACH = 'daint'
    if opts.clm_version == '4.5':
        COMP += 'CLM45'

    out_exe = os.path.abspath(opts.output)
    sourcemods = [os.path.abspath(src_dir) for src_dir in opts.src_mod]

    create_case_fmt = '{:s}/scripts/create_newcase -res {:s} -compset {:s} -mach {:s} -compiler pgi_oas -case {:s}'
    create_case_cmd = create_case_fmt.format(CESM_TRUNK, RES, COMP, MACH, CASEDIR)

    # Build compiling script
    # ----------------------

    with open('build_cesm.bash', mode='w') as script:
        script.write('#!/bin/bash\n')
        script.write('\n')
        script.write('# ----------------------------------------------\n')
        script.write('# Modules\n')
        script.write('# ----------------------------------------------\n')
        script.write('\n')
        if opts.compiler == 'pgi':
            script.write('module switch PrgEnv-cray PrgEnv-pgi\n')
            if opts.compiler_version is not None:
                script.write('module switch pgi pgi/{:s}\n'.format(opts.compiler_version))
        elif opts.compiler == 'intel':
            script.write('module switch PrgEnv-cray PrgEnv-intel\n')
            if opts.compiler_version is not None:
                script.write('module switch intel intel/{:s}\n'.format(opts.compiler_version))
        elif opts.compiler == 'cray' and opts.compiler_version is not None:
            script.write('module switch cce cce/{:s}\n'.format(opts.compiler_version))
        script.write('\n')
        script.write('module load cray-netcdf\n')
        script.write('module load daint-gpu\n')
        script.write('\n')
        script.write('module list\n')
        script.write('\n')
        script.write('# ----------------------------------------------\n')
        script.write('# Create case\n')
        script.write('# ----------------------------------------------\n')
        script.write('\n')
        script.write('{:s}\n'.format(create_case_cmd))
        script.write('\n')
        script.write('# ----------------------------------------------\n')
        script.write('# Setup case\n')
        script.write('# ----------------------------------------------\n')
        script.write('\n')
        script.write('cd {:s}\n'.format(CASEDIR))
        script.write('\n')
        script.write('switch off river routing\n')
        script.write('./xmlchange RTM_MODE="NULL"\n')
        script.write('\n')
        script.write('set transient CO2\n')
        script.write('./xmlchange CCSM_BGC=CO2A,CLM_CO2_TYPE=diagnostic\n')
        if opts.debug:
            script.write('# activate debug mode\n')
            script.write('./xmlchange -file env_build.xml -id DEBUG -val "TRUE"\n')
        script.write('\n')
        script.write('./cesm_setup\n')
        script.write('\n')
        script.write('# ----------------------------------------------\n')
        script.write('# Add source additions/modifications\n')
        script.write('# ----------------------------------------------\n')
        script.write('\n')
        for src_dir in sourcemods:
            print(src_dir)
            for comp in glob('{:s}/src.*'.format(src_dir)):
                print(comp)
                script.write('rsync -avrL {:s} SourceMods\n'.format(comp))
        script.write('\n')
        script.write('# ----------------------------------------------\n')
        script.write('# Build\n')
        script.write('# ----------------------------------------------\n')
        script.write('\n')
        script.write('{:s}.build\n'.format(EXP))
        script.write('rsync -avr bld/cesm.exe {:s}\n'.format(out_exe))

    os.chmod('build_cesm.bash', 0o755)


    # Execute compiling script
    # ------------------------

    if opts.execute:
        check_call(['./build_cesm.bash'])

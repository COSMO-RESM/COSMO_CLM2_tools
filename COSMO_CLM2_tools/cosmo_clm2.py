from __future__ import print_function
from subprocess import check_call
from argparse import ArgumentParser, RawTextHelpFormatter
import f90nml
from datetime import datetime, timedelta
import os
import re
import xml.etree.ElementTree as ET
from glob import glob
from socket import gethostname
import shutil
import warnings

# Date formats
date_fmt_in = '%Y-%m-%d-%H'
date_fmt_cosmo = '%Y%m%d%H'
date_fmt_cesm = '%Y%m%d'


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                           The COSMO-CLM2 case class
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class case(object):
    """Class defining a COSMO-CLM2 case"""

    # Class wide variables
    # ====================
    # Number of tasks per node
    n_tasks_per_node = 12

    # ====
    # Init
    # ====
    def __init__(self, name='COSMO_CLM2', path=None, start_date=None, end_date=None,
                 run_start_date=None, run_length=None, COSMO_exe='./cosmo', CESM_exe='./cesm.exe',
                 submit_script=None, COSMO_archive=None, CESM_archive=None,
                 wall_time='24:00:00', ncosx=None, ncosy=None, ncesm=None, gpu_mode=False,
                 module_purge=False):
        self.name = name
        self.path = path
        self.start_date = start_date
        self.end_date =  end_date
        self.run_start_date = run_start_date
        self.run_length = run_length
        self.COSMO_exe = COSMO_exe
        self.CESM_exe = CESM_exe
        # - ML - Archiving not implemented yet, just a place holder
        self.COSMO_archive = COSMO_archive
        self.CESM_archive = CESM_archive
        self.submit_script = submit_script
        self.wall_time = wall_time
        self.gpu_mode = gpu_mode
        self.module_purge = module_purge
        self._namelists = {'INPUT_ASS': None, 'INPUT_DIA': None, 'INPUT_DYN': None,
                           'INPUT_INI': None, 'INPUT_IO': None, 'INPUT_ORG': None,
                           'INPUT_PHY': None, 'INPUT_SAT': None,
                           'datm_atm_in': None, 'datm_in': None, 'drv_flds_in': None,
                           'drv_in': None, 'lnd_in': None, 'rof_in': None,
                           'atm_modelio.nml': None, 'cpl_modelio.nml': None,
                           'glc_modelio.nml': None, 'ice_modelio.nml': None,
                           'lnd_modelio.nml': None, 'ocn_modelio.nml': None,
                           'rof_modelio.nml': None, 'wav_modelio.nml': None}
        self.ncosx = ncosx
        self.ncosy = ncosy
        self.ncesm = ncesm
        self._n_nodes = None
        self._complete = False

    # Properties
    # ----------
    @property
    def path(self):
        return self._path
    @path.setter
    def path(self, path):
        if path is None:
            self._path = os.path.abspath(os.path.join(os.environ['SCRATCH'], self.name))
        else:
            self._path = os.path.abspath(path)

    @property
    def start_date(self):
        return self._start_date
    @start_date.setter
    def start_date(self, start_date):
        if start_date is not None:
            self._start_date = datetime.strptime(start_date, date_fmt_in)
        else:
            self._start_date = None

    @property
    def end_date(self):
        return self._end_date
    @end_date.setter
    def end_date(self, end_date):
        if end_date is not None:
            self._end_date = datetime.strptime(end_date, date_fmt_in)
        else:
            self._end_date = None

    @property
    def run_start_date(self):
        return self._run_start_date
    @run_start_date.setter
    def run_start_date(self, date):
        if date is not None:
            self._run_start_date = datetime.strptime(date, date_fmt_in)
        else:
            self._run_start_date = self._start_date

    @property
    def run_length(self):
        return self._run_length
    @run_length.setter
    def run_length(self, dt_str):
        self._run_length = dt_str
        # Set run_end_date and runtime
        if dt_str is not None:
            self._run_end_date = min(add_time_from_str(self._run_start_date, dt_str),
                                     self._end_date)
        else:
            self._run_end_date = self._end_date
        self._runtime = self._run_end_date - self._run_start_date

    @property
    def ncosx(self):
        return self._ncosx
    @ncosx.setter
    def ncosx(self, n):
        self._ncosx = n
        if n is not None:
            INPUT_ORG = self.get_namelist('INPUT_ORG')
            INPUT_ORG['runctl']['nprocx'] = n

    @property
    def ncosy(self):
        return self._ncosy
    @ncosy.setter
    def ncosy(self, n):
        self._ncosy = n
        if n is not None:
            INPUT_ORG = self.get_namelist('INPUT_ORG')
            INPUT_ORG['runctl']['nprocy'] = n

    @property
    def ncesm(self):
        return self._ncesm
    @ncesm.setter
    def ncesm(self, n):
        self._ncesm = n
        if n is not None:
            drv_in = self.get_namelist('drv_in')
            for comp in ['atm', 'cpl', 'glc', 'ice', 'lnd', 'ocn', 'rof', 'wav']:
                drv_in['ccsm_pes']['{:s}_ntasks'.format(comp)] = n
            
    # Read-only Properties
    # --------------------
    @property
    def runtime(self):
        return self._runtime

    @property
    def n_nodes(self):
        return self._n_nodes

    @property
    def namelists(self):
        return self._namelists

    @property
    def run_end_date(self):
        return self._run_end_date


    # =======
    # Methods
    # =======
    def get_namelist(self, name):
        if self._namelists[name] is None:
            self._namelists[name] = f90nml.read(os.path.join(self.path, name))
        return self._namelists[name]


    def check_dates(self):
        # - ML - Don't know if this is really usefull
        #        Maybe if initiating case from its path only gets implemented
        INPUT_ORG = self.get_namelist('INPUT_ORG')
        drv_in = self.get_namelist('drv_in')
        start_date_cosmo = datetime.strptime(INPUT_ORG['runctl']['ydate_ini'], date_fmt_cosmo) \
                           + timedelta(hours=INPUT_ORG['runctl']['hstart'])
        start_date_cesm = datetime.strptime(drv_in['seq_timemgr_inparm']['start_ymd'], date_fmt_cesm)
        runtime_cosmo = (INPUT_ORG['runctl']['nstop'] + 1) * INPUT_ORG['runctl']['dt']
        runtime_cesm = drv_in['stop_n']
        
        status = 0
        if start_date_cosmo != start_date_cesm:
            print ("ERROR: start dates are not identical in COSMO and CESM")
            status += 1
        if runtime_cosmo != runtime_cesm:
            print ("ERROR: simulation lengths are not identical in COSMO and CESM")
            status += 1
        if status > 0:
            raise ValueError("start dates and/or simulation lengths do not match between COSMO and CESM")

    
    def create_missing_dirs(self):
        # COSMO
        # -----
        # input
        INPUT_IO = self.get_namelist('INPUT_IO')
        self.mk_miss_path(INPUT_IO['gribin']['ydirini'])
        self.mk_miss_path(INPUT_IO['gribin']['ydirbd'])
        # output
        ngribout = INPUT_IO['ioctl']['ngribout']
        for k in range(ngribout):
            self.mk_miss_path(INPUT_IO['gribout'][k]['ydir'])
        self.mk_miss_path(INPUT_IO['ioctl']['ydir_restart_in'])
        self.mk_miss_path(INPUT_IO['ioctl']['ydir_restart_out'])
        # CESM
        # ----
        # timing
        drv_in = self.get_namelist('drv_in')
        self.mk_miss_path(drv_in['seq_infodata_inparm']['timing_dir'])
        self.mk_miss_path(drv_in['seq_infodata_inparm']['tchkpt_dir'])
        # input
        for name in self._namelists.keys():
            if 'modelio' in name:
                nml = self.get_namelist(name)
                self.mk_miss_path(nml['modelio']['diri'])
                self.mk_miss_path(nml['modelio']['diro'])

                    
    def mk_miss_path(self, rel_path):
        path = os.path.join(self.path, rel_path)
        if not os.path.exists(path):
            os.makedirs(path)
            

    def empty_timing(self):
        drv_in = self.get_namelist('drv_in')
        shutil.rmtree(os.path.join(self.path, drv_in['seq_infodata_inparm']['timing_dir']),
                      ignore_errors=True)
        shutil.rmtree(os.path.join(self.path, drv_in['seq_infodata_inparm']['tchkpt_dir']),
                      ignore_errors=True)
        self.mk_miss_path(drv_in['seq_infodata_inparm']['timing_dir'])
        self.mk_miss_path(drv_in['seq_infodata_inparm']['tchkpt_dir'])
            

    def shift_run_dates(self):
        # Check if simulation is complete
        if self._run_start_date == self._end_date:
            self._complete = True
        else:
            self._run_start_date = self._run_end_date
            # Check if it's the last dummy day for COSMO output
            if self._run_start_date == self._end_date:
                self._run_end_date = self._end_date + timedelta(days=1)
            else:
                self._run_end_date = min(add_time_from_str(self._run_start_date, self._run_length),
                                         self._end_date)
            self._runtime = self._run_end_date - self._run_start_date


    def update_nml_run_dates(self):
        hstart = (self.run_start_date - self.start_date).total_seconds() // 3600.0
        runtime_hours = self.runtime.total_seconds() // 3600.0
        hstop = int(hstart + runtime_hours)
        # INPUT_ORG - runctl
        nml = self.get_namelist('INPUT_ORG')['runctl']
        nml['hstart'] = hstart
        dt = nml['dt']
        nml['nstop'] = int(self.runtime.total_seconds() / dt) - 1
        # INPUT_IO
        gribouts = self.get_namelist('INPUT_IO')['gribout']
        if not isinstance(gribouts, list):
            gribouts = [gribouts]
        for gribout in gribouts:
            gribout['hcomb'][0:2] = hstart, runtime_hours
        self.get_namelist('INPUT_IO')['ioctl']['nhour_restart'] = [hstop, hstop, int(runtime_hours)]
        # drv_in - seq_timemgr_inparm
        nml = self.get_namelist('drv_in')['seq_timemgr_inparm']
        nml['start_ymd'] = int(self.run_start_date.strftime(date_fmt_cesm))
        nml['stop_n'] = int(self.runtime.total_seconds())
        nml['restart_n'] = int(self.runtime.total_seconds())
        if self.run_start_date == self.start_date:
           nml['start_type'] = 'startup'
        else:
           nml['start_type'] = 'continue'
            
    def init_INPUT_ORG(self):
        nml = self.get_namelist('INPUT_ORG')
        nml['runctl']['ydate_ini'] = self.start_date.strftime(date_fmt_cosmo)
        nml['runctl']['hstart'] = (self._run_start_date - self._start_date).total_seconds() / 3600.0
        dt = nml['runctl']['dt']
        nml['runctl']['nstop'] = int(self.runtime.total_seconds() / dt) - 1


    def init_INPUT_IO(self):
        hstart = (self.run_start_date - self.start_date).total_seconds() // 3600.0
        runtime_hours = self.runtime.total_seconds() // 3600.0
        hstop = int(hstart + runtime_hours)
        INPUT_IO = self.get_namelist('INPUT_IO')
        # Output
        if isinstance(INPUT_IO['gribout'], list):
            gribouts_in = INPUT_IO['gribout']
        else:
            gribouts_in = [INPUT_IO['gribout']]
        gribouts_out = []
        for gribout in gribouts_in:
            if runtime_hours >= gribout['hcomb'][2]:
                gribout['hcomb'][0:2] = hstart, runtime_hours
                gribouts_out.append(gribout)
        if gribouts_out:
            INPUT_IO['gribout'] = gribouts_out
            INPUT_IO['ioctl']['ngribout'] = len(gribouts_out)
        else:
            del INPUT_IO['gribout']
        # Restart
        INPUT_IO['ioctl']['nhour_restart'] = [hstop, hstop, int(runtime_hours)]
        

    def update_INPUT_DIA(self):
        # - ML - do domething for M_* files => 1 per run or append?
        raise ValueError('Not implemented yet')

        
    def init_drv_in(self):
        drv_in = self.get_namelist('drv_in')
        drv_in['seq_timemgr_inparm']['start_ymd'] = int(self.run_start_date.strftime(date_fmt_cesm))
        drv_in['seq_timemgr_inparm']['stop_n'] = int(self.runtime.total_seconds())
        drv_in['seq_timemgr_inparm']['restart_n'] = int(self.runtime.total_seconds())
        if self.run_start_date == self.start_date:
           drv_in['seq_timemgr_inparm']['start_type'] = 'startup'
        else:
           drv_in['seq_timemgr_inparm']['start_type'] = 'continue'
        drv_in['seq_infodata_inparm']['case_name'] = self.name
        drv_in['seq_infodata_inparm']['username'] = os.environ['USER']


    def update_namcouple(self):
        with open(os.path.join(self.path, 'namcouple_tmpl'), mode='r') as f:
            content = f.read()
        content = re.sub('_runtime_', str(int(self.runtime.total_seconds())), content)
        with open(os.path.join(self.path, 'namcouple'), mode='w') as f:
            f.write(content)


    def build_proc_config(self):
        nml = self.get_namelist('INPUT_ORG')
        n_cos = nml['runctl']['nprocx'] * nml['runctl']['nprocy']
        nml = self.get_namelist('drv_in')
        n_cesm = nml['ccsm_pes']['atm_ntasks']
        n_tot = n_cos + n_cesm
        # - ML - Add warning if not a round number of nodes
        self._n_nodes = n_tot // self.n_tasks_per_node
        with open(os.path.join(self.path, 'proc_config'), mode='w') as f:
            f.write('{:d}-{:d} ./{:s}\n'.format(0, n_cos-1, os.path.basename(self.COSMO_exe)))
            f.write('{:d}-{:d} ./{:s}\n'.format(n_cos, n_tot-1, os.path.basename(self.CESM_exe)))


    def write_open_nml(self):
        for name, nml in self._namelists.items():
            if nml is not None:
                nml.write(os.path.join(self.path, name), force=True)
                

    def init_submit_script(self):
        account = os.path.normpath(os.environ['PROJECT']).split(os.path.sep)[-2]
        rules = {'#SBATCH +--nodes=.*$': '#SBATCH --nodes={:d}'.format(self.n_nodes),
                 '#SBATCH +--job-name=.*$': '#SBATCH --job-name={:s}'.format(self.name),
                 '#SBATCH +--output=.*$': '#SBATCH --output={:s}.out'.format(self.name),
                 '#SBATCH +--error=.*$': '#SBATCH --error={:s}.out'.format(self.name),
                 '#SBATCH +--account=.*$': '#SBATCH --account={:s}'.format(account),
                 '#SBATCH +--time=.*$': '#SBATCH --time={:s}'.format(self.wall_time)}
        with open(os.path.join(self.path, self.submit_script), mode='r+') as f:
            text = f.read()
            for pattern, repl in rules.items():
                text = re.sub(pattern, repl, text, flags=re.MULTILINE)
            f.seek(0)
            f.write(text)
            f.truncate()


    def build_controller(self):
        account = os.path.normpath(os.environ['PROJECT']).split(os.path.sep)[-2]
        logfile = '{:s}_{:s}-{:s}.out'.format(self.name,
                                              self.run_start_date.strftime(date_fmt_cesm),
                                              self.run_end_date.strftime(date_fmt_cesm))
        with open(os.path.join(self.path, 'controller'), mode='w') as script:
            script.write('#!/bin/bash -l\n')
            script.write('#SBATCH --constraint=gpu\n')
            script.write('#SBATCH --job-name={:s}\n'.format(self.name))
            script.write('#SBATCH --nodes={:d}\n'.format(self.n_nodes))
            script.write('#SBATCH --output={:s}\n'.format(logfile))
            script.write('#SBATCH --error={:s}\n'.format(logfile))
            script.write('#SBATCH --account={:s}\n'.format(account))
            script.write('#SBATCH --time={:s}\n'.format(self.wall_time))
            script.write('\n')
            if self.module_purge:
                script.write('module purge\n')
                script.write('module load PrgEnv-pgi\n')
                script.write('module load cray-netcdf\n')
            else:
                script.write('module switch PrgEnv-cray PrgEnv-pgi\n')
                script.write('module load cray-netcdf\n')
            script.write('module list\n')
            script.write('\n')
            script.write('export MALLOC_MMAP_MAX_=0\n')
            script.write('export MALLOC_TRIM_THRESHOLD_=536870912\n')
            script.write('\n')
            script.write('# Set this to avoid segmentation faults\n')
            script.write('ulimit -s unlimited\n')
            script.write('ulimit -a\n')
            script.write('\n')
            script.write('export OMP_NUM_THREADS=1\n')
            if self.gpu_mode:
                script.write('export MV2_ENABLE_AFFINITY=0\n')
                script.write('export MV2_USE_CUDA=1\n')
                script.write('MPICH_RDMA_ENABLED_CUDA=1\n')
                script.write('export MPICH_G2G_PIPELINE=256\n')
                script.write('\n')
            script.write('cc2_control_case ./config.xml\n')


    def update_controller(self):
        logfile = '{:s}_{:s}-{:s}.out'.format(self.name,
                                              self.run_start_date.strftime(date_fmt_cesm),
                                              self.run_end_date.strftime(date_fmt_cesm))
        rules = {'#SBATCH +--output=.*$': '#SBATCH --output={:s}'.format(logfile),
                 '#SBATCH +--error=.*$': '#SBATCH --error={:s}'.format(logfile)}
        with open(os.path.join(self.path, 'controller'), mode='r+') as f:
            text = f.read()
            for pattern, repl in rules.items():
                text = re.sub(pattern, repl, text, flags=re.MULTILINE)
            f.seek(0)
            f.write(text)
            f.truncate()


    def to_xml(self, file_name):

        def indent(elem, level=0):
            i = "\n" + level*"  "
            if len(elem):
                if not elem.text or not elem.text.strip():
                    elem.text = i + "  "
                if not elem.tail or not elem.tail.strip():
                    elem.tail = i
                for elem in elem:
                    indent(elem, level+1)
                if not elem.tail or not elem.tail.strip():
                    elem.tail = i
            else:
                if level and (not elem.tail or not elem.tail.strip()):
                    elem.tail = i
                    
        config = ET.Element('config')
        tree = ET.ElementTree(config)
        ET.SubElement(config, 'name').text = self.name
        ET.SubElement(config, 'path').text = self.path
        ET.SubElement(config, 'start_date').text = self.start_date.strftime(date_fmt_in)
        ET.SubElement(config, 'end_date').text = self.end_date.strftime(date_fmt_in)
        ET.SubElement(config, 'run_start_date').text = self.run_start_date.strftime(date_fmt_in)
        ET.SubElement(config, 'run_length').text = self.run_length
        ET.SubElement(config, 'COSMO_exe').text = self.COSMO_exe
        ET.SubElement(config, 'CESM_exe').text = self.CESM_exe
        ET.SubElement(config, 'submit_script').text = self.submit_script
        ET.SubElement(config, 'wall_time').text = self.wall_time
        indent(config)
        tree.write(os.path.join(self.path, file_name), xml_declaration=True)


    def submit(self):
        cwd = os.getcwd()
        os.chdir(self.path)
        file_list = glob('YU*') + glob('debug*') + glob('core*')  + glob('nout.*') + glob('*.timers_*')
        for f in file_list:
            os.remove(f)
        check_call(['sbatch', self.submit_script])
        os.chdir(cwd)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                           Module functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def add_time_from_str(date1, dt_str):
    """Increment date from a string

    Return the date resulting from date + N1 years + N2 months or date + N3 days
    where dt_str is a string of the form 'N1yN2m' or 'N1y' or 'N2m' or 'N3d',
    N1, N2 and N3 being arbitrary integers potentially including sign and
    'y', 'm' and 'd' the actual letters standing for year, month and day respectivly."""
        
    ky, km, kd, ny, nm, nd = 0, 0, 0, 0, 0, 0
    for k, c in enumerate(dt_str):
        if c == 'y':
            ky, ny = k, int(dt_str[0:k])
        if c == 'm':
            km, nm = k, int(dt_str[ky:k])
            
    if km == 0 and ky == 0:
        for k, c in enumerate(dt_str):
            if c == 'd':
                kd, nd = k, int(dt_str[0:k])
        if kd == 0:
            raise ValueError("date increment '" + dt_str + "' doesn't have the correct format")
        else:
            return date1 + timedelta(days=nd)
    else:
        y2, m2, d2, h2 = date1.year, date1.month, date1.day, date1.hour
        y2 += ny + (nm+m2-1) // 12
        m2 = (nm+m2-1) % 12 + 1
        return datetime(y2, m2, d2, h2)

        
def case_from_xml(xml_file):
    """Build a COSMO_CLM2 case from xml file"""
    
    config = ET.parse(os.path.normpath(xml_file)).getroot()
    return case(**{opt.tag: opt.text for opt in config.iter() if opt is not config})
        

def create_new_case():
    """Create a new Cosmo-CLM2 case"""

    if "daint" not in gethostname():
        raise ValueError("cosmo_clm2 is only implemented for the Piz Daint machine")

    # Parse setup options from command line and xml file
    # ==================================================
    
    # Options from command line
    # -------------------------
    dsc = "Set up and run a COSMO_CLM2 case\n"\
          "--------------------------------\n"\
          "Options can be set up either by xml file or the following command line arguments.\n"\
          "xml file options must be stored in a subelement of the root element tagged 'cmd_line'.\n"\
          "Command line arguments have precedence over xml file ones."
    parser = ArgumentParser(description=dsc, formatter_class=RawTextHelpFormatter)
    parser.add_argument('--setup-file', metavar='FILE', help="xml file conatining setup options")
    parser.add_argument('--name', help="case name (default: 'COSMO_CLM2')")
    parser.add_argument('--path', help="directory where the case is set up (default: $SCRATCH/NAME)")
    parser.add_argument('--start_date', metavar='DATE_1',
                        help="simulation start date formatted as YYYY-MM-DD-HH")
    parser.add_argument('--end_date', metavar='DATE_2',
                        help="simulation end date formatted as YYYY-MM-DD-HH")
    parser.add_argument('--run_length', metavar='N1yN2m',
                        help="restart every N1 year + N2 month\n"\
                        " N1 and N2 are arbitrary integers potentially including sign\n"\
                        "'y' and 'm' are actual letters standing for 'year' and 'month'\n"\
                        "N1y can be omitted to specify only month (>12 is possible)")
    parser.add_argument('--cos_in', help="COSMO input files directory (default: './COSMO_input')")
    parser.add_argument('--cos_nml', help="COSMO namelists directory (default: './COSMO_nml')")
    parser.add_argument('--cos_exe', help="path to COSMO executable (default: './cosmo')")
    parser.add_argument('--cesm_in', help="CESM input files directory (default: './CESM_input')")
    parser.add_argument('--cesm_nml', help="CESM namelists directory (default: './CESM_nml')")
    parser.add_argument('--cesm_exe', help="CESM executable (default: './cesm.exe')")
    parser.add_argument('--oas_in', help="OASIS input files directory (default: './OASIS_input')")
    parser.add_argument('--oas_nml', help="OASIS namelists directory (default: './OASIS_nml')")
    parser.add_argument('--ncosx', type=int, help="number of subdomains along the 'x-axis'\n"\
                        "for COSMO domain decomposition (default: from INPUT_ORG namelist)")
    parser.add_argument('--ncosy', type=int, help="number of subdomains along the 'y-axis'\n"\
                        "for COSMO domain decomposition (default: from INPUT_ORG namelist)")
    parser.add_argument('--ncesm', type=int, help="number of subdomains for CESM domain decomposition'\n"\
                        " (default: from drv_in namelist)")
    parser.add_argument('--submit_script', help="batch submit script")
    parser.add_argument('--wall_time', help="reserved time on compute nodes (default: '24:00:00')")
    parser.add_argument('--gpu_mode', help="run COSMO on gpu")
    parser.add_argument('--module_purge', help="purge modules before loading and running")
    parser.add_argument('--no_submit', action='store_false', dest='submit',
                        help="Do not submit job after setup\n"\
                        "Only command line argument, cannot be set in xml file")

    opts = parser.parse_args()
    
    # Set options to xml value if needed or default if nothing provided
    # -----------------------------------------------------------------
    defaults = {'name': 'COSMO_CLM2', 'path': None, 'start_date': None, 'end_date': None,
                'run_length': None, 'cos_in': './COSMO_input', 'cos_nml': './COSMO_nml',
                'cos_exe': './cosmo', 'cesm_in': './CESM_input', 'cesm_nml': './CESM_nml',
                'cesm_exe': './cesm.exe', 'oas_in': './OASIS_input', 'oas_nml': './OASIS_nml',
                'submit_script': None, 'wall_time': '24:00:00',
                'ncosx': None, 'ncosy': None, 'ncesm': None, 'gpu_mode': False,
                'module_purge': False}
    if opts.setup_file is not None:
        tree = ET.parse(opts.setup_file)
        xml_node = tree.getroot().find('cmd_line')
    else:
        xml_node = None
    apply_defaults(opts, xml_node, defaults)
    
    # Check if dates are given
    # ------------------------
    if opts.start_date is None or opts.end_date is None:
        raise ValueError("start_date and end_date have to be specified "
                         "either by comand line or xml file")
    
    # Log
    # ===
    log = 'Setting up case {:s} in {:s}'.format(opts.name, opts.path)
    under = '-' * len(log)
    print(log + '\n' + under)
        
    # Transfer data
    # =============
    # - ML - For now, no choice for the i/o directory structure
    if not os.path.exists(opts.path):
        os.makedirs(opts.path)
    transfer_COSMO_input(opts.cos_in, opts.path+'/COSMO_input', opts.start_date, opts.end_date)
    check_call(['rsync', '-avr', opts.cos_nml+'/', opts.path])
    check_call(['rsync', '-avr', opts.cos_exe, opts.path])
    check_call(['rsync', '-avr', opts.cesm_in+'/', opts.path+'/CESM_input/'])
    check_call(['rsync', '-avr', opts.cesm_nml+'/', opts.path])
    check_call(['rsync', '-avr', opts.cesm_exe, opts.path])
    check_call(['rsync', '-avr', opts.oas_in+'/', opts.path])
    check_call(['rsync', '-avr', opts.oas_nml+'/', opts.path])
    if opts.submit_script is not None:
        check_call(['rsync', '-avr', opts.submit_script, opts.path])
        
    # Create case instance
    # ====================
    cc2case = case(name=opts.name, path=opts.path,
                   start_date=opts.start_date, end_date=opts.end_date,
                   run_start_date=opts.start_date, run_length=opts.run_length,
                   COSMO_exe=os.path.basename(opts.cos_exe),
                   CESM_exe=os.path.basename(opts.cesm_exe),
                   submit_script=os.path.basename(opts.submit_script) if opts.submit_script is not None else None,
                   wall_time=opts.wall_time,
                   ncosx=opts.ncosx, ncosy=opts.ncosy, ncesm=opts.ncesm,
                   gpu_mode=opts.gpu_mode,
                   module_purge=opts.module_purge)

    # Init configuration files
    # ========================
    # Change parameters from xml file if required
    if opts.setup_file is not None:
        nodes = tree.getroot().findall('change_par')
        if nodes:
            for node in nodes:
                name = node.get('file')
                block = node.get('block')
                n = node.get('n')
                param = node.get("param")
                value = node.text
                if name is None:
                    raise ValueError("namelist file xml attribute is required to change parameter")
                if block is None:
                    raise ValueError("block xml attribute is required to change parameter")
                if param is None:
                    raise ValueError("param xml attribute is required to change parameter")
                nml = cc2case.get_namelist(name)
                if n is None:
                    nml[block][param] = value
                else:
                    nml[block][int(n)-1][param] = value
                    
    # Init date maters to ensure COSMO and CESM setups match
    cc2case.init_INPUT_ORG()
    cc2case.init_INPUT_IO()
    cc2case.init_drv_in()
    cc2case.update_namcouple()
    cc2case.build_proc_config()

    # Init submit script
    if opts.submit_script is not None:
        cc2case.init_submit_script()
    else:
        cc2case.build_controller()
        
    # Ensure workdir structure exists
    cc2case.create_missing_dirs()

    # Finalize namelists
    cc2case.write_open_nml()

    # Save case configuration to xml file
    # ===================================
    cc2case.to_xml('config.xml')

    # Submit case
    # ===========
    if opts.submit:
        if opts.submit_script is not None:
            cc2case.submit()
        else:
            cwd = os.getcwd()
            os.chdir(cc2case.path)
            check_call(['sbatch', 'controller', './config.xml'])
            os.chdir(cwd)

        
def apply_defaults(opts, xml_node, defaults):
    """Set options with opts > xml_file > defaults"""
    for opt, default  in defaults.items():
        apply_def = False
        if getattr(opts, opt) is None:
            if xml_node is None:
                apply_def = True
            else:
                xml_opt = xml_node.find(opt)
                if xml_opt is None:
                    apply_def = True
                else:
                    opt_val = xml_opt.text
                    if opt_val is None:
                        apply_def = True
                    else:
                        if xml_opt.get('type') is None:
                            setattr(opts, opt, opt_val)
                        else:
                            opt_type = eval(xml_opt.get('type'))
                            if isinstance(opt_type, type):
                                setattr(opts, opt, opt_type(opt_val))
                            else:
                                raise ValueError("xml atribute 'type' for option {:s}".format(opt)
                                                 + " is not a valit python type")
        if apply_def:
            setattr(opts, opt, default)


def transfer_COSMO_input(src_dir, target_dir, start_date, end_date):
    d1 = datetime.strptime(start_date, date_fmt_in)
    # - ML - Add 1 day for the dummy day after simulation end date
    d2 = datetime.strptime(end_date, date_fmt_in) + timedelta(days=1)
    file_list = os.listdir(src_dir)

    with open('transfer_list', mode ='w') as t_list:
        for f in file_list:
            transfer = False
            if f.startswith('laf'):
                transfer = d1 <= datetime.strptime(f[3:13], date_fmt_cosmo) <= d2
            elif f.startswith('lbfd'):
                transfer = d1 <= datetime.strptime(f[4:14], date_fmt_cosmo) <= d2
            if transfer:
                t_list.write(f + '\n')
                
    check_call(['rsync', '-avr', '--files-from', 'transfer_list',
                os.path.normpath(src_dir)+'/', os.path.normpath(target_dir)+'/'])
    
    os.remove('transfer_list')
        

def control_case():
    # Parse arguments
    dsc = "Control a COSMO_CLM2 case"
    parser = ArgumentParser(description=dsc, formatter_class=RawTextHelpFormatter)
    parser.add_argument('xml_path', help="path to xml file containing case description")
    cfg = parser.parse_args()

    # Read case configuration from xml file
    path, xml_file = os.path.split(cfg.xml_path)
    os.chdir(path)
    cc2case = case_from_xml(xml_file)

    # Clean workdir
    file_list = glob('YU*') + glob('debug*') + glob('core*')  + glob('nout.*') + glob('*.timers_*')
    for f in file_list:
        os.remove(f)

    # Run
    check_call(['srun', '-u', '--multi-prog', './proc_config'])

    # Submit next run
    cc2case.shift_run_dates()
    if not cc2case._complete:
        cc2case.update_nml_run_dates()
        cc2case.update_namcouple()
        cc2case.write_open_nml()
        cc2case.empty_timing()
        cc2case.update_controller()
        cc2case.to_xml('config.xml')
        check_call(['sbatch', 'controller', './config.xml'])
        
        
    

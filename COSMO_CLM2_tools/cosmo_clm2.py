# [[file:~/Projects/COSMO_CLM2_tools/COSMO_CLM2_tools.org::*preamble][preamble:1]]
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
import time

# Date formats
date_fmt_in = '%Y-%m-%d-%H'
date_fmt_cosmo = '%Y%m%d%H'
date_fmt_cesm = '%Y%m%d'
# preamble:1 ends here

# [[file:~/Projects/COSMO_CLM2_tools/COSMO_CLM2_tools.org::*class%20case][class case:1]]
class case(object):
    """Class defining a COSMO-CLM2 case"""

    # Class wide variables
    # ====================
    # Number of tasks per node
    _n_tasks_per_node = 12

    # ====
    # Init
    # ====
    def __init__(self, name='COSMO_CLM2', path=None,
                 start_date=None, end_date=None, run_length=None,
                 COSMO_exe='./cosmo', CESM_exe='./cesm.exe',
                 wall_time='24:00:00', account=None,
                 ncosx=None, ncosy=None, ncosio=None, ncesm=None,
                 gpu_mode=False, modules_opt=None, pgi_version=None,
                 dummy_day=True):
        # Basic init (no particular work required)
        self.run_length = run_length
        self.COSMO_exe = COSMO_exe
        self.CESM_exe = CESM_exe
        self.wall_time = wall_time
        self.account = account
        self.gpu_mode = gpu_mode
        self.modules_opt = modules_opt
        self.pgi_version = pgi_version
        self.dummy_day = dummy_day
        # Settings involving namelist changes
        self.path = path
        self.nml = nmldict(self)
        self.name = name
        self.start_date = start_date
        self.end_date = end_date
        self._compute_run_dates()   # defines _run_start_date, _run_end_date and _runtime (maybe _end_date)
        self._apply_run_dates()
        self._check_gribout()
        self.ncosx = ncosx
        self.ncosy = ncosy
        self.ncosio = ncosio
        self.ncesm = ncesm   # Keep that order betwen ncos* and ncesm as ncesm is adapted in gpu mode
        self.write_open_nml()   # Nothing requires changing namelists after that
        # Create batch scripts
        self._build_proc_config()
        self._build_controller()
        # Create missing directories
        self._create_missing_dirs()
        # Write case to xml file
        self.to_xml('config.xml')

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
    def name(self):
        return self._name
    @name.setter
    def name(self, name):
        self._name = name
        self.nml['drv_in']['seq_infodata_inparm']['case_name'] = name

    @property
    def start_date(self):
        return self._start_date
    @start_date.setter
    def start_date(self, start_date):
        if start_date is not None:
            self._start_date = datetime.strptime(start_date, date_fmt_in)
            self.nml['INPUT_ORG']['runctl']['ydate_ini'] = self._start_date.strftime(date_fmt_cosmo)
        elif 'ydate_ini' in self.nml['INPUT_ORG']['runctl'].keys():
            self._start_date = datetime.strptime(self.nml['INPUT_ORG']['runctl']['ydate_ini'],
                                                 date_fmt_cosmo)
        else:
            raise ValueError("ydate_ini has to be given in INPUT_ORG/runctl if no start_date is provided")

    @property
    def end_date(self):
        return self._end_date
    @end_date.setter
    def end_date(self, end_date):
        if end_date is not None:
            self._end_date = datetime.strptime(end_date, date_fmt_in)
            self.nml['INPUT_ORG']['runctl']['ydate_end'] = self._end_date.strftime(date_fmt_cosmo)
        elif 'ydate_end' in self.nml['INPUT_ORG']['runctl'].keys():
            self._end_date = datetime.strptime(self.nml['INPUT_ORG']['runctl']['ydate_end'], date_fmt_cosmo)
        else:
            self._end_date = None

    @property
    def ncosx(self):
        return self._ncosx
    @ncosx.setter
    def ncosx(self, n):
        if n is None:
            self._ncosx = self.nml['INPUT_ORG']['runctl']['nprocx']
        else:
            self._ncosx = n
            self.nml['INPUT_ORG']['runctl']['nprocx'] = n

    @property
    def ncosy(self):
        return self._ncosy
    @ncosy.setter
    def ncosy(self, n):
        if n is None:
            self._ncosy = self.nml['INPUT_ORG']['runctl']['nprocy']
        else:
            self._ncosy = n
            self.nml['INPUT_ORG']['runctl']['nprocy'] = n

    @property
    def ncosio(self):
        return self._ncosio
    @ncosio.setter
    def ncosio(self, n):
        if n is None:
            self._ncosio = self.nml['INPUT_ORG']['runctl']['nprocio']
        else:
            self._ncosio = n
            self.nml['INPUT_ORG']['runctl']['nprocio'] = n

    @property
    def ncesm(self):
        return self._ncesm
    @ncesm.setter
    def ncesm(self, n):
        # total number of COSMO tasks
        self._ncos = self._ncosx * self._ncosy + self._ncosio
        if self.gpu_mode:   # Populate nodes with CESM tasks except one
            self._n_nodes = self._ncos
            self._ncesm = self._n_nodes * (self._n_tasks_per_node - 1)
        else:   # Determine number of CESM tasks and deduce number of nodes
            if n is None:
                self._ncesm = self.nml['drv_in']['ccsm_pes']['lnd_ntasks']
            else:
                self._ncesm = n
            ntot = self._ncos + self._ncesm
            if ntot % self._n_tasks_per_node != 0:
                msg = "total number of tasks (ncosx x ncosy + ncosio + ncesm = {:d}) has to be divisible by {:d}"
                raise ValueError(msg.format(ntot, self._n_tasks_per_node))
            self._n_nodes = ntot // self._n_tasks_per_node
        # Apply number of CESM tasks to all relevant namelist parameters
        for comp in ['atm', 'cpl', 'glc', 'ice', 'lnd', 'ocn', 'rof', 'wav']:
            self.nml['drv_in']['ccsm_pes']['{:s}_ntasks'.format(comp)] = self._ncesm

    @property
    def account(self):
        return self._account
    @account.setter
    def account(self, acc):
        if acc is None:
            # Guess from ${PROJECT} environment variable
            self._account = os.path.normpath(os.environ['PROJECT']).split(os.path.sep)[-2]
        else:
            self._account = acc


    # =======
    # Methods
    # =======
    def _compute_run_dates(self):
        # Access to namelists
        # -------------------
        INPUT_ORG = self.nml['INPUT_ORG']
        drv_in = self.nml['drv_in']
        # Read in _run_start_date
        # -----------------------
        date_cosmo = datetime.strptime(INPUT_ORG['runctl']['ydate_ini'], date_fmt_cosmo) \
                     + timedelta(hours=INPUT_ORG['runctl']['hstart'])
        date_cesm = datetime.strptime(str(drv_in['seq_timemgr_inparm']['start_ymd']), date_fmt_cesm)
        if date_cosmo != date_cesm:
            raise ValueError("start dates are not identical in COSMO and CESM namelists")
        else:
            self._run_start_date = date_cosmo
        # Compute _runtime and _run_end_date (possibly _end_date)
        # -------------------------------------------------------
        if self._end_date is not None:
            if self._run_start_date > self._end_date:
                raise ValueError("run sart date is larger than case end date")
            elif self._run_start_date == self._end_date:
                self._runtime = timedelta(days=1)
                self._run_end_date = self._end_date + self._runtime
            else:
                if self.run_length is None:
                    self._run_end_date = self._end_date
                else:
                    self._run_end_date = min(add_time_from_str(self._run_start_date, self.run_length),
                                             self._end_date)
                self._runtime = self._run_end_date - self._run_start_date
        else:
            if self.run_length is None:
                runtime_cosmo = (INPUT_ORG['runctl']['nstop'] + 1) * INPUT_ORG['runctl']['dt'] \
                                - INPUT_ORG['runctl']['hstart'] * 3600.0
                runtime_cesm = drv_in['seq_timemgr_inparm']['stop_n']
                if runtime_cosmo != runtime_cesm:
                    raise ValueError("run lengths are not identical in COSMO and CESM namelists")
                else:
                    self._runtime = timedelta(seconds=runtime_cosmo)
                    self._run_end_date = self._run_start_date + self._runtime
            else:
                self._run_end_date = add_time_from_str(self._run_start_date, self.run_length)
                self._runtime = self._run_end_date - self._run_start_date
            self._end_date = self._run_end_date


    def _apply_run_dates(self):
        # Compute times
        hstart = (self._run_start_date - self.start_date).total_seconds() // 3600.0
        runtime_seconds = self._runtime.total_seconds()
        runtime_hours = runtime_seconds // 3600.0
        hstop = hstart + runtime_hours
        # Access to namelists
        INPUT_ORG = self.nml['INPUT_ORG']
        INPUT_IO = self.nml['INPUT_IO']
        drv_in = self.nml['drv_in']
        # adapt INPUT_ORG
        INPUT_ORG['runctl']['nstop'] = int(hstop * 3600.0 // INPUT_ORG['runctl']['dt']) - 1
        # adapt INPUT_IO
        for gribout in self._get_gribouts():
            gribout['hcomb'][0:2] = hstart, hstop
        INPUT_IO['ioctl']['nhour_restart'] = [int(hstop), int(hstop), 24]
        # adapt drv_in
        drv_in['seq_timemgr_inparm']['stop_n'] = int(runtime_seconds)
        drv_in['seq_timemgr_inparm']['restart_n'] = int(runtime_seconds)
        # adapt namcouple
        with open(os.path.join(self.path, 'namcouple_tmpl'), mode='r') as f:
            content = f.read()
        content = re.sub('_runtime_', str(int(self._runtime.total_seconds())), content)
        with open(os.path.join(self.path, 'namcouple'), mode='w') as f:
            f.write(content)


    def _check_gribout(self):
        # Only keep gribout blocks that fit within runtime
        # (essentially to avoid crash for short tests)
        runtime_hours = self._runtime.total_seconds() // 3600.0
        gribouts_out = []
        gribouts_in = self._get_gribouts()
        for gribout in gribouts_in:
            if runtime_hours >= gribout['hcomb'][2]:
                gribouts_out.append(gribout)
        if gribouts_out:
            self.nml['INPUT_IO']['gribout'] = gribouts_out
            self.nml['INPUT_IO']['ioctl']['ngribout'] = len(gribouts_out)
        else:
            if gribouts_in:
                del self.nml['INPUT_IO']['gribout']


    def _get_gribouts(self):
        if 'gribout' not in self.nml['INPUT_IO'].keys():
            return []
        else:
            gribouts = self.nml['INPUT_IO']['gribout']
            if not isinstance(gribouts, list):
                gribouts = [gribouts]
            return gribouts


    def write_open_nml(self):
        self.nml.write_all()


    def _create_missing_dirs(self):
        # COSMO
        # -----
        # input
        self._mk_miss_path(self.nml['INPUT_IO']['gribin']['ydirini'])
        self._mk_miss_path(self.nml['INPUT_IO']['gribin']['ydirbd'])
        # output
        for gribout in self._get_gribouts():
            self._mk_miss_path(gribout['ydir'])
        self._mk_miss_path(self.nml['INPUT_IO']['ioctl']['ydir_restart_in'])
        self._mk_miss_path(self.nml['INPUT_IO']['ioctl']['ydir_restart_out'])
        # CESM
        # ----
        # timing
        # - ML - remove if exists before creating
        shutil.rmtree(os.path.join(self.path, self.nml['drv_in']['seq_infodata_inparm']['timing_dir']),
                      ignore_errors=True)
        shutil.rmtree(os.path.join(self.path, self.nml['drv_in']['seq_infodata_inparm']['tchkpt_dir']),
                      ignore_errors=True)
        self._mk_miss_path(self.nml['drv_in']['seq_infodata_inparm']['timing_dir'])
        self._mk_miss_path(self.nml['drv_in']['seq_infodata_inparm']['tchkpt_dir'])
        # input / output
        for comp in ['atm', 'cpl', 'glc', 'ice', 'lnd', 'ocn', 'rof', 'wav']:
            self._mk_miss_path(self.nml['{:s}_modelio.nml'.format(comp)]['modelio']['diri'])
            self._mk_miss_path(self.nml['{:s}_modelio.nml'.format(comp)]['modelio']['diro'])


    def _mk_miss_path(self, rel_path):
        path = os.path.join(self.path, rel_path)
        if not os.path.exists(path):
            print('Creating path ' + path)
            os.makedirs(path)


    def _build_proc_config(self):
        with open(os.path.join(self.path, 'proc_config'), mode='w') as f:
            if self.gpu_mode:
                N = self._n_tasks_per_node
                line = ",".join([str(k*N) for k in range(self._n_nodes)])
                f.write("{:s} ./{:s}\n".format(line, self.COSMO_exe))
                line = ",".join(["{:d}-{:d}".format(k*N+1,(k+1)*N-1) for k in range(self._n_nodes)])
                f.write("{:s} ./{:s}\n".format(line, self.CESM_exe))
            else:
                f.write('{:d}-{:d} ./{:s}\n'.format(0, self._ncos-1, self.COSMO_exe))
                f.write('{:d}-{:d} ./{:s}\n'.format(self._ncos, self._ncos+self._ncesm-1, self.CESM_exe))


    def _build_controller(self):
        logfile = '{:s}_{:s}-{:s}.out'.format(self.name,
                                              self._run_start_date.strftime(date_fmt_cesm),
                                              self._run_end_date.strftime(date_fmt_cesm))
        with open(os.path.join(self.path, 'controller'), mode='w') as script:
            script.write('#!/usr/bin/env bash\n')
            script.write('#SBATCH --constraint=gpu\n')
            script.write('#SBATCH --job-name={:s}\n'.format(self.name))
            script.write('#SBATCH --nodes={:d}\n'.format(self._n_nodes))
            script.write('#SBATCH --output={:s}\n'.format(logfile))
            script.write('#SBATCH --error={:s}\n'.format(logfile))
            script.write('#SBATCH --account={:s}\n'.format(self.account))
            script.write('#SBATCH --time={:s}\n'.format(self.wall_time))
            script.write('\n')
            script.write('module list\n')
            if self.gpu_mode:
                script.write('source /etc/bash.bashrc.local\n')
            if self.modules_opt:
                if self.modules_opt == 'purge':
                    script.write('module purge\n')
                    script.write('module load PrgEnv-pgi\n')
                elif self.modules_opt == 'switch':
                    script.write('module switch PrgEnv-cray PrgEnv-pgi\n')
                if self.pgi_version == '16.9.0':
                    script.write('module unload pgi\n')
                    script.write('module load pgi/16.9.0\n')
                elif self.pgi_version == '17.10.0':
                    script.write('module unload pgi\n')
                    script.write('module use /apps/common/UES/pgi/17.10/modulefiles\n')
                    script.write('module load pgi/17.10\n')
                    script.write('export PGI_VERS_STR=17.10.0\n')
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
                script.write('\n')
                script.write('# Use for gpu mode\n')
                script.write('export MV2_ENABLE_AFFINITY=0\n')
                script.write('export MV2_USE_CUDA=1\n')
                script.write('MPICH_RDMA_ENABLED_CUDA=1\n')
                script.write('export MPICH_G2G_PIPELINE=256\n')
                script.write('\n')
            script.write('cc2_control_case ./config.xml\n')


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
        ET.SubElement(config, 'run_length').text = self.run_length
        ET.SubElement(config, 'COSMO_exe').text = self.COSMO_exe
        ET.SubElement(config, 'CESM_exe').text = self.CESM_exe
        ET.SubElement(config, 'wall_time').text = self.wall_time
        ET.SubElement(config, 'account').text = self.account
        ET.SubElement(config, 'gpu_mode', attrib={'type': 'bool'}).text = '1' if self.gpu_mode else ''
        if self.modules_opt is not None:
            ET.SubElement(config, 'modules_opt').text = self.modules_opt
        if self.pgi_version is not None:
            ET.SubElement(config, 'pgi_version').text = self.pgi_version
        ET.SubElement(config, 'dummy_day', attrib={'type': 'bool'}).text = '1' if self.dummy_day else ''
        indent(config)
        tree.write(os.path.join(self.path, file_name), xml_declaration=True)


    def set_next_run(self):
        if ((self._run_start_date >= self._end_date) or
            (self._run_end_date == self._end_date and not self.dummy_day)):
            return False
        else:
            hstart = (self._run_end_date - self._start_date).total_seconds() // 3600.0
            self.nml['INPUT_ORG']['runctl']['hstart'] = hstart
            self.nml['drv_in']['seq_timemgr_inparm']['start_ymd'] = int(self._run_end_date.strftime(date_fmt_cesm))
            self._compute_run_dates()
            # - ML - Setting ydirini might not be needed, try without at some point
            self.nml['INPUT_IO']['gribin']['ydirini'] = self.nml['INPUT_IO']['ioctl']['ydir_restart_out']
            for gribout in self._get_gribouts():
                gribout['lwrite_const'] = False
            self.nml['drv_in']['seq_infodata_inparm']['start_type'] = 'continue'
            self.write_open_nml()
            self._update_controller()
            return True


    def _update_controller(self):
        logfile = '{:s}_{:s}-{:s}.out'.format(self.name,
                                              self._run_start_date.strftime(date_fmt_cesm),
                                              self._run_end_date.strftime(date_fmt_cesm))
        rules = {'#SBATCH +--output=.*$': '#SBATCH --output={:s}'.format(logfile),
                 '#SBATCH +--error=.*$': '#SBATCH --error={:s}'.format(logfile)}
        with open(os.path.join(self.path, 'controller'), mode='r+') as f:
            content = f.read()
            for pattern, repl in rules.items():
                content = re.sub(pattern, repl, content, flags=re.MULTILINE)
            f.seek(0)
            f.write(content)
            f.truncate()


    def submit(self):
        cwd = os.getcwd()
        os.chdir(self.path)
        check_call(['sbatch', 'controller', './config.xml'])
        os.chdir(cwd)


    def run(self):
        cwd = os.getcwd()
        # Clean workdir
        os.chdir(self.path)
        file_list = glob('YU*') + glob('debug*') + glob('core*')  + glob('nout.*') + glob('*.timers_*')
        for f in file_list:
            os.remove(f)
        # Run
        start_time = time.time()
        check_call(['srun', '-u', '--multi-prog', './proc_config'])
        elapsed = time.time() - start_time
        print("\nCase {name:s} ran in {elapsed:.2f}\n".format(name=self.name, elapsed=elapsed))
        os.chdir(cwd)
# class case:1 ends here

# [[file:~/Projects/COSMO_CLM2_tools/COSMO_CLM2_tools.org::*class%20nmldict][class nmldict:1]]
class nmldict(dict):
    """Dictionnary of all the namelists of a case. Only load tha namelist if needed"""
    def __init__(self, cc2case):
        dict.__init__(self)
        self.cc2case = cc2case

    def __getitem__(self, key):
        if key not in self:
            self[key] = f90nml.read(os.path.join(self.cc2case.path, key))
        return dict.__getitem__(self, key)

    def write(self, name):
        self[name].write(os.path.join(self.cc2case.path, name), force=True)

    def write_all(self):
        for name, nml in self.items():
            self.write(name)
# class nmldict:1 ends here

# [[file:~/Projects/COSMO_CLM2_tools/COSMO_CLM2_tools.org::*add_time_from_str][add_time_from_str:1]]
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
# add_time_from_str:1 ends here

# [[file:~/Projects/COSMO_CLM2_tools/COSMO_CLM2_tools.org::*case_from_xml][case_from_xml:1]]
def case_from_xml(xml_file):
    """Build a COSMO_CLM2 case from xml file"""

    config = ET.parse(os.path.normpath(xml_file)).getroot()
    args={}
    for opt in config.iter():
        if opt is not config:
            if opt.get('type') is None:
                args[opt.tag] = opt.text
            else:
                opt_type = eval(opt.get('type'))
                if isinstance(opt_type, type):
                    args[opt.tag] = opt_type(opt.text)
                else:
                    raise ValueError("xml atribute 'type' for option {:s}".format(opt.tag)
                                     + " is not a valid python type")

    return case(**args)
# case_from_xml:1 ends here

# [[file:~/Projects/COSMO_CLM2_tools/COSMO_CLM2_tools.org::*create_new_case][create_new_case:1]]
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
    parser.add_argument('-s', '--setup-file', metavar='FILE', help="xml file conatining setup options")
    parser.add_argument('--name', help="case name (default: 'COSMO_CLM2')")
    parser.add_argument('--path', help="directory where the case is set up (default: $SCRATCH/NAME)")
    parser.add_argument('--start_date', metavar='DATE_1',
                        help="simulation start date formatted as YYYY-MM-DD-HH")
    parser.add_argument('--end_date', metavar='DATE_2',
                        help="simulation end date formatted as YYYY-MM-DD-HH")
    parser.add_argument('--run_length', metavar='dt',
                        help="sets simulation length if end_date not specified or run length\n"\
                        "between restarts otherwise\n"\
                        "dt is of the form 'N1yN2m' or 'N1y' or 'N2m' or 'N3d'\n"\
                        "N1, N2 and N3 being arbitrary integers (N2>12 possible) and\n"\
                        "'y', 'm' and 'd' stand for year, month and day")
    parser.add_argument('--cos_in', help="COSMO input files directory (default: './COSMO_input')")
    parser.add_argument('--cos_nml', help="COSMO namelists directory (default: './COSMO_nml')")
    parser.add_argument('--cos_exe', help="path to COSMO executable (default: './cosmo')")
    parser.add_argument('--cesm_in', help="CESM input files directory (default: './CESM_input')")
    parser.add_argument('--cesm_nml', help="CESM namelists directory (default: './CESM_nml')")
    parser.add_argument('--cesm_exe', help="CESM executable (default: './cesm.exe')")
    parser.add_argument('--oas_in', help="OASIS input files directory (default: './OASIS_input')")
    parser.add_argument('--oas_nml', help="OASIS namelists directory (default: './OASIS_nml')")
    parser.add_argument('--ncosx', type=int, help="number of subdomains along the 'x-axis'\n"\
                        "for COSMO domain decomposition (type: int, default: from INPUT_ORG namelist)")
    parser.add_argument('--ncosy', type=int, help="number of subdomains along the 'y-axis'\n"\
                        "for COSMO domain decomposition (type: int, default: from INPUT_ORG namelist)")
    parser.add_argument('--ncosio', type=int, help="number of cores dedicated to i/o work'\n"\
                        "(type: int, default: from INPUT_ORG namelist)")
    parser.add_argument('--ncesm', type=int, help="number of subdomains for CESM domain decomposition'\n"\
                        "(type: int, default: from drv_in namelist)")
    parser.add_argument('--wall_time', help="reserved time on compute nodes (default: '24:00:00')")
    parser.add_argument('--account', help="account to use for batch script (default: infered from $PROJECT)")
    parser.add_argument('--gpu_mode', type=bool, help="run COSMO on gpu (type: bool, default: False)")
    parser.add_argument('--modules_opt', choices=['switch', 'purge'],
                        help="Option for loading modules at run time (default: None)")
    parser.add_argument('--pgi_version', choices=['16.9.0', '17.5.0', '17.10.0'],
                        help="specify pgi compiler version (default: None)")
    parser.add_argument('--dummy_day', type=bool,
                        help="perform a dummy day run after end of simulation to get last COSMO output.\n"\
                        "(type: bool, default: True)")
    parser.add_argument('--no_submit', action='store_false', dest='submit',
                        help="do not submit job after setup\n"\
                        "only command line argument, cannot be set in xml file")
    parser.add_argument('--gen_oasis', action='store_true',
                        help="generate OASIS auxiliary files\n"\
                        "note that OASIS will crash after producing the files\n"\
                        "only command line argument, cannot be set in xml file\n"
                        )

    opts = parser.parse_args()
    if opts.gen_oasis:
        opts.dummy_day = False

    # Set options to xml value if needed or default if nothing provided
    # -----------------------------------------------------------------
    defaults = {'name': 'COSMO_CLM2', 'path': None, 'start_date': None, 'end_date': None,
                'run_length': None, 'cos_in': './COSMO_input', 'cos_nml': './COSMO_nml',
                'cos_exe': './cosmo', 'cesm_in': './CESM_input', 'cesm_nml': './CESM_nml',
                'cesm_exe': './cesm.exe', 'oas_in': './OASIS_input', 'oas_nml': './OASIS_nml',
                'ncosx': None, 'ncosy': None, 'ncosio': None, 'ncesm': None,
                'wall_time': '24:00:00', 'account': None, 'dummy_day': True,
                'gpu_mode': False, 'modules_opt': None, 'pgi_version': None}
    if opts.setup_file is not None:
        tree = ET.parse(opts.setup_file)
        xml_node = tree.getroot().find('cmd_line')
    else:
        xml_node = None
    apply_defaults(opts, xml_node, defaults)

    if opts.path is None:
        opts.path = os.path.join(os.environ['SCRATCH'], opts.name)

    # Log
    # ===
    log = 'Setting up case {:s} in {:s}'.format(opts.name, opts.path)
    under = '-' * len(log)
    print(log + '\n' + under)

    # Transfer data
    # =============
    # - ML - For now, no choice for the I/O directory structure
    # - ML - Do first transfering namelists, then create case, then transfer input
    if not os.path.exists(opts.path):
        os.makedirs(opts.path)
    INPUT_IO = f90nml.read(os.path.join(opts.cos_nml, 'INPUT_IO'))
    dh = INPUT_IO['gribin']['hincbound']
    ext =''
    if 'yform_read' in INPUT_IO['ioctl'].keys():
        if INPUT_IO['ioctl']['yform_read'] == 'ncdf':
            ext = '.nc'
    transfer_COSMO_input(opts.cos_in, opts.path+'/COSMO_input',
                         opts.start_date, opts.end_date,
                         opts.run_length, dh, opts.dummy_day, ext)
    check_call(['rsync', '-avr', opts.cos_nml+'/', opts.path])
    check_call(['rsync', '-avr', opts.cos_exe, opts.path])
    check_call(['rsync', '-avr', opts.cesm_in+'/', opts.path+'/CESM_input/'])
    check_call(['rsync', '-avr', opts.cesm_nml+'/', opts.path])
    check_call(['rsync', '-avr', opts.cesm_exe, opts.path])
    if not opts.gen_oasis:
        check_call(['rsync', '-avr', opts.oas_in+'/', opts.path])
    else:
        for f in os.listdir(opts.oas_in):
            os.remove(os.path.join(opts.path, f))
    check_call(['rsync', '-avr', opts.oas_nml+'/', opts.path])

    # Create case instance
    # ====================
    cc2case = case(name=opts.name, path=opts.path,
                   start_date=opts.start_date, end_date=opts.end_date,
                   run_length=opts.run_length,
                   COSMO_exe=os.path.basename(opts.cos_exe),
                   CESM_exe=os.path.basename(opts.cesm_exe),
                   wall_time=opts.wall_time, account=opts.account,
                   ncosx=opts.ncosx, ncosy=opts.ncosy, ncosio=opts.ncosio, ncesm=opts.ncesm,
                   gpu_mode=opts.gpu_mode,
                   modules_opt=opts.modules_opt, pgi_version=opts.pgi_version,
                   dummy_day=opts.dummy_day)

    # Change parameters from xml file if required
    # ===========================================
    # Change namelist parameters from xml file
    if opts.setup_file is not None:
        nodes = tree.getroot().findall('change_par')
        if nodes:
            for node in nodes:
                name = node.get('file')
                block = node.get('block')
                n = node.get('n')
                param = node.get("param")
                val_str = node.text
                if name is None:
                    raise ValueError("namelist file xml attribute is required to change parameter")
                if block is None:
                    raise ValueError("block xml attribute is required to change parameter")
                if param is None:
                    raise ValueError("param xml attribute is required to change parameter")
                nml = cc2case.nml[name][block]
                if node.get('type') is None:
                    value = val_str
                elif node.get('type') == 'py_eval':
                    value = eval(val_str)
                else:
                    val_type = eval(node.get('type'))
                    if isinstance(val_type, type):
                        value = val_type(val_str)
                    else:
                        err_mess = "Given xml atribute 'type' for parameter {:s} is {:s}\n"\
                                   "It has to be either 'py_eval' or a valid build in python type"
                        raise ValueError(err_mess.format(param, val_type))
                if n is None:
                    nml[param] = value
                else:
                    nml[int(n)-1][param] = value
    # Change namelist parameters from certain cmd line arguments
    if opts.gen_oasis:
        cc2case.nml['drv_in']['ccsm_pes']['atm_ntasks'] = 1

    # Finalize
    # ========
    cc2case.write_open_nml()
    cc2case.to_xml('config.xml')

    # Submit case
    # ===========
    if opts.submit:
        cc2case.submit()
# create_new_case:1 ends here

# [[file:~/Projects/COSMO_CLM2_tools/COSMO_CLM2_tools.org::*apply_defaults][apply_defaults:1]]
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
                    opt_val_str = xml_opt.text
                    if xml_opt.get('type') is None:
                        setattr(opts, opt, opt_val_str)
                    elif xml_opt.get('type') == 'py_eval':
                        setattr(opts, opt, eval(opt_val_str))
                    else:
                        opt_type = eval(xml_opt.get('type'))
                        if isinstance(opt_type, type):
                            setattr(opts, opt, opt_type(opt_val_str))
                        else:
                            raise ValueError("xml atribute 'type' for option {:s}".format(opt)
                                             + " is not a valid python type")
        if apply_def:
            setattr(opts, opt, default)
# apply_defaults:1 ends here

# [[file:~/Projects/COSMO_CLM2_tools/COSMO_CLM2_tools.org::*transfer_COSMO_input][transfer_COSMO_input:1]]
def transfer_COSMO_input(src_dir, target_dir, start_date, end_date,
                         run_length, dh, dummy_day, ext):

    d1 = datetime.strptime(start_date, date_fmt_in)
    if end_date is None:
        if run_length is None:
            raise ValueError("if end_date is none, provide run_length")
        else:
            d2 = add_time_from_str(d1, run_length)
    else:
        d2 = datetime.strptime(end_date, date_fmt_in)
    delta = timedelta(seconds=dh*3600.0)

    def check_input(root, date, file_list, dummy=False):
        file_name = root + format(date.strftime(date_fmt_cosmo)) + ext
        if os.path.exists(os.path.join(src_dir, file_name)):
            file_list.write(file_name + '\n')
            return True
        elif dummy:
            raise ValueError("Creating dummy input files: no tool available on Piz Daint\n"\
                             "to alter the date of the input file, wether grib or netcdf.\n"\
                             "Please proceed manually.")
            # dummy_date = datetime(d1.year, d1.month, d1.day, date.hour)
            # dummy_file_name = root + format(dummy_date.strftime(date_fmt_cosmo)) + ext
            # msg = "WARNING: Copying {:s} as {:s} for additionnal dummy day (produce last COSMO output)"
            # print(msg.format(dummy_file_name, file_name))
            # in_file = os.path.join(src_dir, dummy_file_name)
            # out_file = os.path.join(target_dir, file_name)
            # if ext == '':
            #     check_call(['grib_set', '-s', 'dataDate={:s}'.format(dummy_date.strftime('%Y%m%d')),
            #                 in_file, out_file])
            # else:
            #     raise ValueError("Creating dummy input files: no tool available on Piz Daint\n"\
            #                      "to alter the date of netcdf file, please proceed manually.")
            #     # shutil.copy(in_file, out_file)
            # return False
        else:
            raise ValueError("input file {:s} is missing".format(file_name))

    # Check all input files for current period
    with open('transfer_list', mode ='w') as t_list:
        check_input('laf', d1, t_list)
        cur_date = d1
        while cur_date <= d2:
            check_input('lbfd', cur_date, t_list)
            cur_date += delta
    check_call(['rsync', '-avr', '--files-from', 'transfer_list',
                os.path.normpath(src_dir)+'/', os.path.normpath(target_dir)+'/'])

    # Add a dummy day to produce last COSMO output
    if dummy_day:
        do_transfer = False
        with open('transfer_list', mode ='w') as t_list:
            while cur_date <= d2 + timedelta(days=1):
                do_transfer = do_transfer or check_input('lbfd', cur_date, t_list, dummy=True)
                cur_date += delta
        if do_transfer:
            check_call(['rsync', '-avr', '--files-from', 'transfer_list',
                        os.path.normpath(src_dir)+'/', os.path.normpath(target_dir)+'/'])

    os.remove('transfer_list')
# transfer_COSMO_input:1 ends here

# [[file:~/Projects/COSMO_CLM2_tools/COSMO_CLM2_tools.org::*control_case][control_case:1]]
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

    # Run
    cc2case.run()

    # Submit next run
    if cc2case.set_next_run():
        cc2case.submit()
# control_case:1 ends here

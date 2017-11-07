from subprocess import check_call
from argparse import ArgumentParser, RawTextHelpFormatter
import f90nml
from datetime import datetime, timedelta
import os
import re
import xml.etree.ElementTree as ET
import fileinput
from glob import glob
from socket import gethostname
import shutil

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                           The case class
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class case():
    """Class defining a COSMO-CLM2 case"""

    # Class wide variables
    # ====================
    # Date formats
    date_fmt_in = '%Y-%m-%d-%H'
    date_fmt_cosmo = '%Y%m%d%H'
    date_fmt_cesm = '%Y%m%d'

    # Number of tasks per node
    n_tasks_per_node = 12

    # ====
    # Init
    # ====
    def __init__(self, name='COSMO_CLM2', path=None, start_date=None, end_date=None,
                 run_start_date=None, run_length=None,
                 COSMO_input='./COSMO_input', CESM_input='./CESM_input',
                 COSMO_exe='./cosmo', CESM_exe='./cesm.exe', submit_script='./submit',
                 COSMO_archive=None, CESM_archive=None):
        self.name = name
        self.path = path
        self.start_date = start_date
        self.end_date =  end_date
        self.run_start_date = run_start_date
        self.run_length = run_length
        self.COSMO_input = COSMO_input
        self.CESM_input = CESM_input
        self.COSMO_exe = COSMO_exe
        self.CESM_exe = CESM_exe
        self.COSMO_archive=COSMO_archive
        self.CESM_archive=CESM_archive
        self.submit_script = submit_script
        self._namelists = {'COSMO': {'INPUT_ASS': None, 'INPUT_DIA': None, 'INPUT_DYN': None,
                                     'INPUT_INI': None, 'INPUT_IO': None, 'INPUT_ORG': None,
                                     'INPUT_PHY': None, 'INPUT_SAT': None},
                           'CESM': {'datm_atm_in': None, 'datm_in': None, 'drv_flds_in': None,
                                    'drv_in': None, 'lnd_in': None, 'rof_in': None,
                                    'atm_modelio.nml': None, 'cpl_modelio.nml': None,
                                    'glc_modelio.nml': None, 'ice_modelio.nml': None,
                                    'lnd_modelio.nml': None, 'ocn_modelio.nml': None,
                                    'rof_modelio.nml': None, 'wav_modelio.nml': None}}
        self._n_nodes = None

    # Properties
    # ----------
    @property
    def path(self):
        return self._path

    @path.setter
    def path(self, path):
        if path is None:
            self._path = os.path.abs_path(os.path.join(os.environ['SCRATCH'], self.name))
        else:
            self._path = os.path.abs_path(path)

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
    def run_start_date(self, run_start_date):
        if run_start_date is not None:
            self._run_start_date = datetime.strptime(run_start_date, date_fmt_in)
        else:
            self._run_start_date = self._start_date

    @property
    def run_length(self):
        return self._run_length

    @run_length.setter
    def run_length(self, dt_str):
        if dt_str is not None:
            self._run_end_date = min(add_time_from_str(self._run_start_date, dt_str),
                                     self.end_date)
        else:
            self._run_end_date = self._end_date
        self._runtime = self._run_end_date - self._run_start_date
            
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
    def get_namelist(self, model, name):
        if self._namelists[model][name] is None:
            self._namelists[model][name] = f90nml.read(os.path.join(self.path, name))
        return = self._namelists[model][name]


    def check_dates(self):
        # - ML - Don't know if this is really usefull
        #        Maybe if initiating case from its path only gets implemented
        INPUT_ORG = self.get_namelist('COSMO', 'INPUT_ORG')
        drv_in = self.get_namelist('CESM', 'drv_in')
        date_ini_cosmo = 
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

    
    def create_workdir_structure(self):
        os.makedirs(os.path.join(self.path, 'COSMO_output'), exist_ok=True)
        os.makedirs(os.path.join(self.path, 'CESM_timing'), exist_ok=True)


    def update_nml_run_dates(self):
        nml = self.get_namelist('COSMO', 'INPUT_ORG')
        nml['runctl']['hstart'] = (self.run_start_date - self.start_date).total_seconds() / 3600.0
        dt = nml['runctl']['dt']
        nml['runctl']['nstop'] = int(self.runtime.total_seconds() / dt) - 1
        nml = self.get_namelist('CESM', 'drv_in')
        nml['seq_timemgr_inparm']['start_ymd'] = int(self.run_start_date.strftime(date_fmt_cesm))
        nml['seq_timemgr_inparm']['stop_n'] = int(self.runtime.total_seconds())
        nml['seq_timemgr_inparm']['restart_n'] = int(self.runtime.total_seconds())


    def set_next_run_dates(self):
        self._run_start_date = self._run_end_date
        self._run_end_date = add_time_from_str(self._run_start_date, self._run_length)
        self._runtime = self._run_end_date - self._run_start_date
        
        
    def update_INPUT_ORG(self):
        nml = self.get_namelist('COSMO', 'INPUT_ORG')
        nml['runctl']['ydate_ini'] = self.start_date.strftime(date_fmt_cosmo)
        nml['runctl']['hstart'] = (self._run_start_date - self._start_date).total_seconds() / 3600.0
        dt = nml['runctl']['dt']
        nml['runctl']['nstop'] = int(self.runtime.total_seconds() / dt) - 1


    def update_INPUT_IO(self):
        nml = self.get_namelist('COSMO', 'INPUT_IO')
        runtime_hours = self.runtime.total_seconds() // 3600.0
        # Input
        nml['gribin']['ydirini'] = self.COSMO_input
        nml['gribin']['ydirbd'] = self.COSMO_input
        # Output
        if isinstance(nml['gribout'], list):
            gribouts_in = nml['gribout']
        else:
            gribouts_in = [nml['gribout']]
        gribouts_out = []
        for gribout in gribouts_in:
            if runtime_hours >= gribout['hcomb'][2]:
                gribout['hcomb'][1] = float(runtime_hours)
                gribout['ydir'] = './COSMO_output/'
                gribouts_out.append(gribout)
        if gribouts_out:
            nml['gribout'] = gribouts_out
        else:
            del nml['gribout']
        # Restart
        nml['ioctl']['nhour_restart'] = [runtime_hours, runtime_hours, runtime_hours]
        nml['ioctl']['ydir_restart_in'] = './COSMO_output/'
        nml['ioctl']['ydir_restart_out'] = './COSMO_output/'
        

    def update_INPUT_DIA(self):
        # - ML - do domething for M_* files => 1 per run or append?

        
    def update_drv_in(self):
        nml = self.get_namelist('CESM', 'drv_in')
        nml['seq_timemgr_inparm']['start_ymd'] = int(self.run_start_date.strftime(date_fmt_cesm))
        nml['seq_timemgr_inparm']['stop_n'] = int(self.runtime.total_seconds())
        nml['seq_timemgr_inparm']['restart_n'] = int(self.runtime.total_seconds())
        nml['seq_infodata_inparm']['case_name'] = CASE_NAME
        nml['seq_infodata_inparm']['username'] = os.environ['USER']

        
    def update_modelio(self):
        for name in self._namelists['CESM'].keys():
            if 'modelio' in name:
                nml = self.get_namelist('CESM', name)
                nml['modelio']['diri'] = self.CESM_input
                nml['modelio']['diro'] = '.'


    def update_namcouple(self):
        shutil.copyfile('namcouple_tmpl', 'namcouple')
        with fileinput.FileInput(os.path.join(self.path, 'namcouple'), inplace=True) as namcouple:
            for line in namcouple:
                print(line.replace('_runtime_', str(int(self.runtime.total_seconds()))), end='')


    def build_proc_config(self):
        nml = self.get_namelist('COSMO', 'INPUT_ORG')
        n_cos = nml['runctl']['nprocx'] * nml['runctl']['nprocy']
        nml = self.get_namelist('CESM', 'drv_in')
        n_cesm = nml['ccsm_pes']['atm_ntasks']
        n_tot = n_cos + n_cesm
        # - ML - Add warning if not a round number of nodes
        self._n_nodes = n_tot // n_tasks_per_node
        with open(os.path.join(self.path, 'proc_config'), mode='w') as f:
            f.write('{:d}-{:d} ./{:s}\n'.format(0, n_cos-1, os.path.basename(self.COSMO_exe)))
            f.write('{:d}-{:d} ./{:s}\n'.format(n_cos, n_tot-1, os.path.basename(self.CESM_exe)))


    def update_submit_script(self):
        account = os.path.normpath(os.environ['PROJECT']).split(os.path.sep)[-2]
        rules = {'#SBATCH +--nodes=.*$': '#SBATCH --nodes={:d}'.format(self.n_nodes),
                 '#SBATCH +--job-name=.*$': '#SBATCH --job-name={:s}'.format(self.name),
                 '#SBATCH +--output=.*$': '#SBATCH --output={:s}.out'.format(self.name),
                 '#SBATCH +--account=.*$': '#SBATCH --account={:s}.out'.format(account)}
        with open(self.submit_script, mode='r+') as f:
            text = f.read()
            for pattern, repl in rules.items():
                text = re.sub(pattern, repl, text, flags=re.MULTILINE)
            f.seek(0)
            f.write(text)
            f.truncate()


    def write_open_nml(self):
        for model in _namelists.keys():
            for name, nml in _namelists[model].items():
                if nml is not None:
                    nml.write(os.path.join(self.path, name), force=True)


    def to_xml(self):

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
        ET.SubElement(config, 'run_length').text = self.run_run_length
        ET.SubElement(config, 'COSMO_input').text = self.COSMO_input
        ET.SubElement(config, 'CESM_input').text = self.CESM_input
        ET.SubElement(config, 'COSMO_exe').text = self.COSMO_exe
        ET.SubElement(config, 'CESM_exe').text = self.CESM_exe
        ET.SubElement(config, 'submit_script').text = self.submit_script
        indent(config)
        tree.write(os.path.join(self.path, 'config.xml'), xml_declaration=True)


    def submit(self, clean=True):
        cwd = os.getcwd()
        os.chdir(self.path)
        if clean:
            file_list = glob('YU*') + glob(self.name+'*')
            for f in file_list:
                os.remove(f)
        check_call(['sbatch', self.submit_script])
        os.chdir(cwd)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                           Module functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def add_time_from_str(date, dt_str):
    """Return the date resulting from date + n1 years + n2 months
    where dt_str has the form 'n1yn2m', n1 and n2 being arbitrary integers
    potentially including sign and 'y' and 'm' the actual letters
    standing for year and month respectivly"""
    
    ky, ny, nm = 0, 0 , 0
    for k, c in enumerate(dt_str):
        if char == 'y':
            ky, ny = k, int(dt_str[0:k])
            if char == 'm':
                nm = int(dt_str[ky:k])
                date = self._run_start_date
        y1, m1, d1, h1 = date.year, date.month, date.day, date.hour
        y2 = y1 + ny + (nm+m1-1) // 12
        m2 = (nm+m1-1) % 12
        return datetime(y2, m2, d1, h1)

        
def case_from_xml(xml_file):
    """Build a COSMO_CLM2 case from xml file"""
    
    tree = ET.parse(os.path.normpath(xml_file))
    config = tree.getroot()
    options = {}
    for opt in config.iter():
        options[opt.tag] = opt.text

    return case(**options)
        

def create_new_case():
    """Create a new Cosmo-CLM2 case"""

    if "daint" not in gethostname():
        raise ValueError("cosmo_clm2 is only implemented for the Piz Daint machine")

    # Parse setup options from command line and xml file
    # ==================================================
    
    # Options from command line
    # -------------------------
    dsc = "Set up a COSMO_CLM2 case.\n"\
          "\n"\
          "Options can be set up either by xml file\n"\
          "or the following command line arguments"
    parser = ArgumentParser(description=dsc, formatter_class=RawTextHelpFormatter)
    parser.add_argument('--setup-file', metavar='FILE', default=None,
                        help="xml file conatining setup options")
    parser.add_argument('--name', default='COSMO_CLM2',
                        help="case name (default: 'COSMO_CLM2')")
    parser.add_argument('--path', default=None,
                        help="directory where the case is set up (default: $SCRATCH/NAME)")
    parser.add_argument('--start_date', metavar='DATE_1',
                        help="simulation start date formatted as YYYY-MM-DD-HH")
    parser.add_argument('--end_date', , metavar='DATE_2',
                        help="simulation end date formatted as YYYY-MM-DD-HH")
    parser.add_argument('--restart_every', metavar='N1yN2m'
                        help="run length: N1 and N2 are arbitrary integers potentially including sign,\n"\
                        "'y' and 'm' are actual letters standing for 'year' and 'month'")
    parser.add_argument('--cos_in', default='./COSMO_input',
                        help="COSMO input files directory (default: './COSMO_input')")
    parser.add_argument('--cesm_in', default='./CESM_input',
                        help="CESM input files directory (default: './CESM_input')")
    parser.add_argument('--oas_in', default='./OASIS_input',
                        help="OASIS input files directory (default: './OASIS_input')")
    parser.add_argument('--cos_nml',default='./COSMO_nml',
                        help="COSMO namelists directory (default: './COSMO_nml')")
    parser.add_argument('--cesm_nml',default='./CESM_nml',
                        help="CESM namelists directory (default: './CESM_nml')")
    parser.add_argument('--oas_nml',default='./OASIS_nml',
                        help="OASIS namelists directory (default: './OASIS_nml')")
    parser.add_argument('--cos_exe',default='./cosmo',
                        help="path to COSMO executable (default: './cosmo')")
    parser.add_argument('--cesm_exe',default='./cesm.exe',
                        help="CESM executable (default: './cesm.exe')")
    parser.add_argument('--submit_script',default='./submit',
                        help="batch submit script (default: './submit')")
    parser.add_argument('--no_submit', action='store_false',
                        help="Do not submit job after setup\n"\
                        "Only command line argument, cannot be set in xml file")

    cmd_opts = vars(parser.parse_args())

    # Options from xml file
    # ---------------------
    xml_opts = {}
    tree = ET.parse(os.path.normpath(xml_file))
    cmd_line = tree.getroot().find('cmd_line')
    xml_opts = {opt.tag: opt.text for opt in cmd_line.iter() if opt is not cmd_line}

    # Final option list
    # -----------------
    # Command line take precedence over xml file
    setup_opts = xml_opts
    for opt, val in cmd_opts.items():
        setup_opts[opt] = val
    if 'start_date' not in setup_opts or 'end_date' not in setup_opts:
        raise ValueError("start_date and end_date have to be specified "
                         "either by comand line or xml file")
    if setup_opts['path'] is None:
        setup_opts['path'] = os.path.join(os.environ['SCRATCH'], setup_opts['name'])
    if 'run_length' not in setup_opts:
        setup_opts['run_length'] = None

    # Transfer data
    # =============
    # - ML - For now, no choice for the i/o directory structure
    CASE_path = setup_opts[path]
    os.makedirs(CASE_path, exist_ok=True)
    check_call(['rsync', '-avr', setup_opts['cos_in']+'/', CASE_path+'/COSMO_input/'])
    check_call(['rsync', '-avr', setup_opts['cos_nml']+'/', CASE_path])
    check_call(['rsync', '-avr', setup_opts['cos_exe'], CASE_path])
    check_call(['rsync', '-avr', setup_opts['cesm_in']+'/', CASE_path+'/CESM_input/'])
    check_call(['rsync', '-avr', setup_opts['cesm_nml']+'/', CASE_path])
    check_call(['rsync', '-avr', setup_opts['cesm_exe'], CASE_path])
    check_call(['rsync', '-avr', setup_opts['oas_in']+'/', CASE_path])
    check_call(['rsync', '-avr', setup_opts['oas_nml']+'/', CASE_path])
    check_call(['rsync', '-avr', setup_opts['submit_script'], CASE_PATH])
        
    # Create case instance
    # ====================
    case_opts = {'name': setup_opts['name'], 'path': setup_opts[path],
                 'start_date': setup_opts['start_date'], 'end_date': setup_opts['end_date'],
                 'run_start_date': setup_opts['start_date'],
                 'COSMO_exe': os.path.basename(setup_opts['cos_exe']),
                 'CESM_exe': os.path.basename(setup_opts['cesm_exe']),
                 'submit_script': os.path.basename(setup_opts['submit_script'])}
    if 'restart_every' in setup_opts:
        case_opts['run_length'] = setup_opts['restart_every']
    cc2case = case(**case_opts)

    # Update files
    # ============
    cc2case.create_workdir_structure()
    cc2case.update_INPUT_ORG()
    cc2case.update_INPUT_IO()
    cc2case.update_drv_in()
    cc2case.update_modelio()
    cc2case.update_namcouple()
    cc2case.build_proc_config()
    cc2case.update_submit_script()
    cc2case.write_open_nml()

    # Save case configuration to xml file
    # ===================================
    cc2case.to_xml()

    # Submit case
    # ===========
    if not setup_opts['no_submit']:
        cc2case.submit()

from __future__ import print_function
from .date_formats import date_fmt_in, date_fmt_cosmo, date_fmt_cesm
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
import sys
from warnings import warn

class base_case(object):
    """Base class defining a COSMO-CLM2 case"""

    _n_tasks_per_node = None
    NotImplementMessage = "required method {:s} not implemented by class {:s}.\n" \
                          "Implement with a single pass statement if irrelevant to this machine."


    def __init__(self, name='COSMO_CLM2', path=None,
                 start_date=None, end_date=None, run_length=None,
                 COSMO_exe='./cosmo', CESM_exe='./cesm.exe',
                 ncosx=None, ncosy=None, ncosio=None, ncesm=None,
                 gpu_mode=False, dummy_day=True, cosmo_only=False,
                 gen_oasis=False):
        # Basic init (no particular work required)
        self.cosmo_only = cosmo_only
        self.gen_oasis = gen_oasis
        self.run_length = run_length
        self.COSMO_exe = COSMO_exe
        if not self.cosmo_only:
            self.CESM_exe = CESM_exe
        self.gpu_mode = gpu_mode
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
        self._organize_tasks(ncosx, ncosy, ncosio, ncesm)
        self.write_open_nml()   # Nothing requires changing namelists after that
        # Create batch script
        self._build_controller()
        # Create missing directories
        self._create_missing_dirs()

    @property
    def path(self):
        return self._path
    @path.setter
    def path(self, path):
        self._path = os.path.abspath(path)

    @property
    def name(self):
        return self._name
    @name.setter
    def name(self, name):
        self._name = name
        if not self.cosmo_only:
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
    def account(self):
        return self._account
    @account.setter
    def account(self, acc):
        if acc is None:
            # Guess from ${PROJECT} environment variable
            self._account = os.path.normpath(os.environ['PROJECT']).split(os.path.sep)[-2]
        else:
            self._account = acc


    def _organize_tasks(self, ncosx, ncosy, ncosio, ncesm):
        # COSMO tasks
        # -----------
        if ncosx is None:
            self._ncosx = self.nml['INPUT_ORG']['runctl']['nprocx']
        else:
            self._ncosx = ncosx
            self.nml['INPUT_ORG']['runctl']['nprocx'] = ncosx
        if ncosy is None:
            self._ncosy = self.nml['INPUT_ORG']['runctl']['nprocy']
        else:
            self._ncosy = ncosy
            self.nml['INPUT_ORG']['runctl']['nprocy'] = ncosy
        if ncosio is None:
            self._ncosio = self.nml['INPUT_ORG']['runctl']['nprocio']
        else:
            self._ncosio = ncosio
            self.nml['INPUT_ORG']['runctl']['nprocio'] = ncosio
        self._ncos = self._ncosx * self._ncosy + self._ncosio
        # CESM tasks and number of nodes
        # ------------------------------
        if self.cosmo_only:
            self._ncesm = 0
            if self.gpu_mode:
                self._n_nodes = self._ncos
            else:
                self._n_nodes = self._ncos // self._n_tasks_per_node
        else:
            if self.gpu_mode:   # Populate nodes with CESM tasks except one
                self._n_nodes = self._ncos
                self._ncesm = self._n_nodes * (self._n_tasks_per_node - 1)
            else:   # Determine number of CESM tasks and deduce number of nodes
                if ncesm is None:
                    self._ncesm = self.nml['drv_in']['ccsm_pes']['lnd_ntasks']
                else:
                    self._ncesm = ncesm
                ntot = self._ncos + self._ncesm
                if ntot % self._n_tasks_per_node != 0:
                    msg = "total number of tasks (ncosx x ncosy + ncosio + ncesm = {:d}) has to be divisible by {:d}"
                    raise ValueError(msg.format(ntot, self._n_tasks_per_node))
                self._n_nodes = ntot // self._n_tasks_per_node
            # Apply number of CESM tasks to all relevant namelist parameters
            for comp in ['atm', 'cpl', 'glc', 'ice', 'lnd', 'ocn', 'rof', 'wav']:
                self.nml['drv_in']['ccsm_pes']['{:s}_ntasks'.format(comp)] = self._ncesm
            if self.gen_oasis:
                self.nml['drv_in']['ccsm_pes']['atm_ntasks'] = 1


    def _compute_run_dates(self):
        # Access to namelists
        # -------------------
        INPUT_ORG = self.nml['INPUT_ORG']
        if not self.cosmo_only:
            drv_in = self.nml['drv_in']
        # Read in _run_start_date
        # -----------------------
        date_cosmo = datetime.strptime(INPUT_ORG['runctl']['ydate_ini'], date_fmt_cosmo) \
                     + timedelta(hours=INPUT_ORG['runctl']['hstart'])
        if not self.cosmo_only:
            date_cesm = datetime.strptime(str(drv_in['seq_timemgr_inparm']['start_ymd']), date_fmt_cesm)
            if date_cosmo != date_cesm:
                raise ValueError("start dates are not identical in COSMO and CESM namelists")
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
                if not self.cosmo_only:
                    runtime_cesm = drv_in['seq_timemgr_inparm']['stop_n']
                    if runtime_cosmo != runtime_cesm:
                        raise ValueError("run lengths are not identical in COSMO and CESM namelists")
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
        if not self.cosmo_only:
            drv_in = self.nml['drv_in']
        # adapt INPUT_ORG
        INPUT_ORG['runctl']['nstop'] = int(hstop * 3600.0 // INPUT_ORG['runctl']['dt']) - 1
        # adapt INPUT_IO
        for gribout in self._get_gribouts():
            gribout['hcomb'][0:2] = hstart, hstop
        INPUT_IO['ioctl']['nhour_restart'] = [int(hstop), int(hstop), 24]
        if not self.cosmo_only:
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
        if not self.cosmo_only:
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


    def _build_controller(self):
        """Place holder for _build_controller method to be implemented by machine specific classes."""

        raise NotImplementedError(NotImplementMessage.format('_build_controller(self)', self.__class__.__name__))


    def _update_xml_config(self, config):
        """Place holder for _update_xml_config method to be implemented by machine specific classes."""

        raise NotImplementedError(NotImplementMessage.format('_update_xml_config(self)', self.__class__.__name__))


    def to_xml(self, file_name='config.xml'):

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
        ET.SubElement(config, 'cosmo_only', type='py_eval').text = str(self.cosmo_only)
        ET.SubElement(config, 'gen_oasis', type='py_eval').text = str(self.gen_oasis)
        ET.SubElement(config, 'start_date').text = self.start_date.strftime(date_fmt_in)
        ET.SubElement(config, 'end_date').text = self.end_date.strftime(date_fmt_in)
        ET.SubElement(config, 'run_length').text = self.run_length
        ET.SubElement(config, 'COSMO_exe').text = self.COSMO_exe
        if not self.cosmo_only:
            ET.SubElement(config, 'CESM_exe').text = self.CESM_exe
        ET.SubElement(config, 'gpu_mode', type='py_eval').text = str(self.gpu_mode)
        ET.SubElement(config, 'dummy_day', type='py_eval').text = str(self.dummy_day)
        self._update_xml_config(config)
        indent(config)
        tree.write(os.path.join(self.path, file_name), xml_declaration=True)


    @staticmethod
    def from_xml(xml_file):
        """Build a case from xml file"""

        config = ET.parse(os.path.normpath(xml_file)).getroot()
        machine_node = config.find('machine')
        if machine_node is None:
            raise ValueError("machine node not found in {:s}".format(xml_file))
        else:
            machine = machine_node.text

        args={}
        for opt in config.iter():
            if opt is not config and opt is not machine_node:
                if opt.get('type') is None:
                    args[opt.tag] = opt.text
                elif opt.get('type') == 'py_eval':
                    args[opt.tag] = eval(opt.text)
                else:
                    opt_type = eval(opt.get('type'))
                    if isinstance(opt_type, type):
                        args[opt.tag] = opt_type(opt.text)
                    else:
                        raise ValueError("xml atribute 'type' for option {:s}".format(opt.tag)
                                         + " has to be a valid python or 'py_eval'")

        if machine == 'daint':
            return daint_case(**args)
        else:
            raise NotImplementedError("machine {:s} not implemeted".format(machine))


    def set_next_run(self):
        if ((self._run_start_date >= self._end_date) or
            (self._run_end_date == self._end_date and not self.dummy_day)):
            return False
        else:
            hstart = (self._run_end_date - self._start_date).total_seconds() // 3600.0
            self.nml['INPUT_ORG']['runctl']['hstart'] = hstart
            if not self.cosmo_only:
                self.nml['drv_in']['seq_timemgr_inparm']['start_ymd'] = int(self._run_end_date.strftime(date_fmt_cesm))
            self._compute_run_dates()
            # - ML - Setting ydirini might not be needed, try without at some point
            self.nml['INPUT_IO']['gribin']['ydirini'] = self.nml['INPUT_IO']['ioctl']['ydir_restart_out']
            for gribout in self._get_gribouts():
                gribout['lwrite_const'] = False
            if not self.cosmo_only:
                self.nml['drv_in']['seq_infodata_inparm']['start_type'] = 'continue'
            self.write_open_nml()
            self._update_controller()
            return True


    def _update_controller(self):
        """Place holder for _update_controller method to be implemented by machine specific classes."""

        raise NotImplementedError(NotImplementMessage.format('_update_controller(self)', self.__class__.__name__))


    def submit(self):
        cwd = os.getcwd()
        os.chdir(self.path)
        self._submit_func()
        os.chdir(cwd)


    def run(self):
        start_time = time.time()
        cwd = os.getcwd()

        # Clean workdir
        os.chdir(self.path)
        file_list = glob('YU*') + glob('debug*') + glob('core*') + glob('nout.*') + glob('*.timers_*')
        for f in file_list:
            os.remove(f)

        # Run
        self._run_func()

        os.chdir(cwd)
        elapsed = time.time() - start_time
        print("\nCase {name:s} ran in {elapsed:.2f}\n".format(name=self.name, elapsed=elapsed))


    def _run_func(self):
        """Place holder for _run_func method to be implemented by machine specific classes."""

        raise NotImplementedError(NotImplementMessage.format('_run_func(self)', self.__class__.__name__))


    def _submit_func(self):
        """Place holder for _submit_func method to be implemented by machine specific classes."""

        raise NotImplementedError(NotImplementMessage.format('_submit_func(self)', self.__class__.__name__))

class nmldict(dict):
    """Dictionnary of all the namelists of a case. Only load the namelist if needed"""
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

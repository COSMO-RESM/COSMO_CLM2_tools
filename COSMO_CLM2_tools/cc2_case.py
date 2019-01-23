from __future__ import print_function
from .tools import date_fmt, add_time_from_str, COSMO_input_file_name
from subprocess import check_call
from argparse import ArgumentParser, RawTextHelpFormatter
import f90nml
from datetime import datetime, timedelta
import os
import re
import xml.etree.ElementTree as ET
from glob import glob
import shutil
import time
import sys

available_cases = {}

def factory(machine, **case_args):
    if machine not in available_cases:
        raise ValueError("machine {:s} not available".format(machine))
    else:
        return available_cases[machine](**case_args)

class cc2_case(object):
    """Base class defining a COSMO-CLM2 case"""

    _target_machine = None
    _n_tasks_per_node = None
    _default_install_dir = None
    NotImplementMessage = "required method {:s} not implemented by class {:s}.\n" \
                          "Implement with a single pass statement if irrelevant to this machine."


    def __init__(self, name='COSMO_CLM2', install_dir=None, install=False,
                 cos_nml='./COSMO_nml', cos_in='./COSMO_input', cos_exe='./cosmo', cos_in_file_size=None,
                 cesm_nml='./CESM_nml', cesm_in='./CESM_input', cesm_exe='./cesm.exe',
                 oas_in='./OASIS_input', oas_nml='./OASIS_nml',
                 start_date=None, end_date=None, run_length=None,
                 ncosx=None, ncosy=None, ncosio=None, ncesm=None,
                 gpu_mode=False, dummy_day=True, cosmo_only=False,
                 gen_oasis=False, input_type='file'):

        # Basic init (no particular work required)
        self.run_length = run_length
        self.gpu_mode = gpu_mode
        self.dummy_day = dummy_day
        self.cosmo_only = cosmo_only
        self.gen_oasis = gen_oasis
        self.cos_in_file_size = cos_in_file_size
        self.cos_in = os.path.abspath(cos_in)
        self.install = install
        self.input_type = input_type
        # Create namelists dictionnary
        self.nml = nmldict(self)
        # Set case name, install_dir and path
        self._name = name
        self.install_dir = install_dir   # also sets self._path
        # Install: transfer namelists, executables and input files
        if self.install:
            log = 'Setting up case {:s} in {:s}'.format(self._name, self._path)
            print(log + '\n' + '-' * len(log))
            self.install_case(cos_nml, cos_in, cos_exe, cesm_nml, cesm_in, cesm_exe, oas_nml, oas_in)
            # Setting case name in namelist not possible before actually transfering the namelists
            self.nml['drv_in']['seq_infodata_inparm']['case_name'] = self.name
        self.cos_exe = cos_exe
        if not self.cosmo_only:
            self.cesm_exe = cesm_exe
        # Settings involving namelist changes
        self.start_date = start_date
        self.end_date = end_date
        self._compute_run_dates()   # defines _run_start_date, _run_end_date and _runtime (maybe _end_date)
        self._apply_run_dates()
        self._check_gribout()
        self._organize_tasks(ncosx, ncosy, ncosio, ncesm)
        # Finish install
        if self.install:
            # Transfer COSMO input files
            # - ML - Change this in future versions: only transfer the first chunck of input files
            self.transfer_cos_in(self.start_date, self.end_date)
            # Create batch script
            # - ML - for later: probably no need to store executable names in config.xml
            self._build_controller()
            # Create missing directories
            self._create_missing_dirs()
        # Check presence and size of input files
        # - ML - for later: only for current chunck
        self._check_COSMO_input(self.start_date, self.end_date)
        # Write namelists object to file
        self.write_open_nml()

    @property
    def cos_exe(self):
        return self._cos_exe
    @cos_exe.setter
    def cos_exe(self, exe_path):
        self._cos_exe = os.path.basename(exe_path)

    @property
    def cesm_exe(self):
        return self._cesm_exe
    @cesm_exe.setter
    def cesm_exe(self, exe_path):
        self._cesm_exe = os.path.basename(exe_path)

    @property
    def install_dir(self):
        return self._install_dir
    @install_dir.setter
    def install_dir(self, ins_dir):
        if ins_dir is None:
            if self._default_install_dir is None:
                raise NotImplementedError("_default_install_dir class variable not set for Class {:s}".format(cls.__name__))
            else:
                self._install_dir = self._default_install_dir
        else:
            self._install_dir = ins_dir
        # Make install_dir absolute
        self._install_dir = os.path.abspath(self._install_dir)
        # Set case path
        self._path = os.path.join(self._install_dir, self.name)

    @property
    def path(self):
        return self._path

    @property
    def name(self):
        return self._name

    @property
    def start_date(self):
        return self._start_date
    @start_date.setter
    def start_date(self, start_date):
        if start_date is not None:
            self._start_date = datetime.strptime(start_date, date_fmt['in'])
            self.nml['INPUT_ORG']['runctl']['ydate_ini'] = self._start_date.strftime(date_fmt['cosmo'])
        elif 'ydate_ini' in self.nml['INPUT_ORG']['runctl'].keys():
            self._start_date = datetime.strptime(self.nml['INPUT_ORG']['runctl']['ydate_ini'],
                                                 date_fmt['cosmo'])
        else:
            raise ValueError("ydate_ini has to be given in INPUT_ORG/runctl if no start_date is provided")

    @property
    def end_date(self):
        return self._end_date
    @end_date.setter
    def end_date(self, end_date):
        if end_date is not None:
            self._end_date = datetime.strptime(end_date, date_fmt['in'])
            self.nml['INPUT_ORG']['runctl']['ydate_end'] = self._end_date.strftime(date_fmt['cosmo'])
        elif 'ydate_end' in self.nml['INPUT_ORG']['runctl'].keys():
            self._end_date = datetime.strptime(self.nml['INPUT_ORG']['runctl']['ydate_end'], date_fmt['cosmo'])
        else:
            self._end_date = None


    def install_case(self, cos_nml, cos_in, cos_exe, cesm_nml, cesm_in, cesm_exe, oas_nml, oas_in):
        if not os.path.exists(self.path):
            # Create case directory
            os.makedirs(self.path)

        # Transfer everything except COSMO input files
        check_call(['rsync', '-avrL', os.path.abspath(cos_nml)+'/', self.path])
        check_call(['rsync', '-avrL', os.path.abspath(cos_exe), self.path])
        if not self.cosmo_only:
            if self.input_type == 'symlink':
                check_call(['ln', '-sf', os.path.abspath(cesm_in), os.path.join(self.path,'CESM_input')])
            elif self.input_type == 'file':
                check_call(['rsync', '-avrL', os.path.abspath(cesm_in)+'/', os.path.join(self.path,'CESM_input')+'/'])
            check_call(['rsync', '-avrL', os.path.abspath(cesm_nml)+'/', self.path])
            check_call(['rsync', '-avrL', os.path.abspath(cesm_exe), self.path])
            if not self.gen_oasis:
                if self.input_type == 'symlink':
                    check_call(['ln', '-sf', os.path.abspath(oas_in), self.path])
                elif self.input_type == 'file':
                    check_call(['rsync', '-avrL', os.path.abspath(oas_in)+'/', self.path])
            else:
                print('generate OASIS file:')
                for f in os.listdir(oas_in):
                    try:
                        print('   removing ' +  os.path.join(self.path, f))
                        os.remove(os.path.join(self.path, f))
                    except OSError:
                        pass
            check_call(['rsync', '-avrL', os.path.abspath(oas_nml)+'/', self.path])


    def _cos_input_delta_ext(self):
        # Set time interval between 2 intput files
        delta = timedelta(hours=self.nml['INPUT_IO']['gribin']['hincbound'])
        # Set file extension
        ext = ''
        if 'yform_read' in self.nml['INPUT_IO']['ioctl']:
            if self.nml['INPUT_IO']['ioctl']['yform_read'] == 'ncdf':
                ext = '.nc'
        return delta, ext


    def transfer_cos_in(self, start_date, end_date):
        delta, ext = self._cos_input_delta_ext()
        if self.install:
            file_name = COSMO_input_file_name('lbfd', self.start_date, ext)
            file_path = os.path.join(self.cos_in, file_name)
            self.cos_in_file_size = os.stat(file_path).st_size
            if self.input_type == 'symlink':
                try:
                    os.makedirs(os.path.join(self.path,'COSMO_input'))
                except OSError:
                    pass
        # function to check and add file to transfer list or directly symlink
        def check_add_file(root, date, file_list):
            file_name = COSMO_input_file_name(root, date, ext)
            if os.path.exists(os.path.join(self.cos_in, file_name)):
                if self.input_type == 'symlink':
                    check_call(['ln', '-sf', os.path.join(self.cos_in, file_name),
                                os.path.join(self.path,'COSMO_input')])
                elif self.input_type == 'file':
                    file_list.write(file_name + '\n')
            else:
                raise ValueError("input file {:s} is missing from {:s}".format(file_name, self.cos_in))
        # Build file list to transfer or symlink
        with open('transfer_list', mode ='w') as t_list:
            check_add_file('laf', start_date, t_list)
            cur_date = start_date
            while cur_date <= end_date:
                check_add_file('lbfd', cur_date, t_list)
                cur_date += delta
        # Tranfer files
        if self.input_type == 'file':
            check_call(['rsync', '-avrL', '--files-from', 'transfer_list',
                        self.cos_in+'/', os.path.join(self.path,'COSMO_input')+'/'])
        # Remove transfer list
        os.remove('transfer_list')


    def _check_COSMO_input(self, start_date, end_date):
        delta, ext = self._cos_input_delta_ext()
        cur_date = start_date
        while cur_date <= end_date:
            file_name = COSMO_input_file_name('lbfd', cur_date, ext)
            file_path = os.path.join(self.path, 'COSMO_input', file_name)
            if not os.path.exists(file_path):
                raise ValueError("COSMO input file {:s} missing".format(file_name))
            fs = os.stat(file_path).st_size
            if fs != self.cos_in_file_size:
                err_mess = "COSMO input file {:s} has byte size {:d} instead of {:d}"
                raise ValueError(err_mess.format(file_name, fs, self.cos_in_file_size))
            cur_date += delta


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
        date_cosmo = datetime.strptime(INPUT_ORG['runctl']['ydate_ini'], date_fmt['cosmo']) \
                     + timedelta(hours=INPUT_ORG['runctl']['hstart'])
        if not self.cosmo_only:
            date_cesm = datetime.strptime(str(drv_in['seq_timemgr_inparm']['start_ymd']), date_fmt['cesm'])
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
        if 'hstop' in INPUT_ORG['runctl']:
            del INPUT_ORG['runctl']['hstop']
        INPUT_ORG['runctl']['nstop'] = int(hstop * 3600.0 // INPUT_ORG['runctl']['dt']) - 1
        if 'hstop' in INPUT_ORG['runctl']:
            del INPUT_ORG['runctl']['hstop']
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
        if 'ydir_restart_in' in self.nml['INPUT_IO']['ioctl']:
            self._mk_miss_path(self.nml['INPUT_IO']['ioctl']['ydir_restart_in'])
        if 'ydir_restart_out' in self.nml['INPUT_IO']['ioctl']:
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

        config_node = ET.Element('config')
        tree = ET.ElementTree(config_node)
        ET.SubElement(config_node, 'machine').text = self._target_machine
        main_node = ET.SubElement(config_node, 'main')
        machine_node = ET.SubElement(config_node, self._target_machine)
        ET.SubElement(main_node, 'name').text = self.name
        ET.SubElement(main_node, 'install_dir').text = self.install_dir
        ET.SubElement(main_node, 'cosmo_only', type='py_eval').text = str(self.cosmo_only)
        ET.SubElement(main_node, 'gen_oasis', type='py_eval').text = str(self.gen_oasis)
        ET.SubElement(main_node, 'start_date').text = self.start_date.strftime(date_fmt['in'])
        ET.SubElement(main_node, 'end_date').text = self.end_date.strftime(date_fmt['in'])
        ET.SubElement(main_node, 'run_length').text = self.run_length
        ET.SubElement(main_node, 'cos_exe').text = self.cos_exe
        ET.SubElement(main_node, 'cos_in_file_size', type='int').text = str(self.cos_in_file_size)
        if not self.cosmo_only:
            ET.SubElement(main_node, 'cesm_exe').text = self.cesm_exe
        ET.SubElement(main_node, 'cos_in').text = self.cos_in
        ET.SubElement(main_node, 'gpu_mode', type='py_eval').text = str(self.gpu_mode)
        ET.SubElement(main_node, 'dummy_day', type='py_eval').text = str(self.dummy_day)
        self._update_xml_config(machine_node)
        indent(config_node)
        tree.write(os.path.join(self.path, file_name), xml_declaration=True)


    def set_next_run(self):
        if ((self._run_start_date >= self._end_date) or
            (self._run_end_date == self._end_date and not self.dummy_day)):
            continue_run = False
        else:
            continue_run = True
            hstart = (self._run_end_date - self._start_date).total_seconds() // 3600.0
            self.nml['INPUT_ORG']['runctl']['hstart'] = hstart
            if not self.cosmo_only:
                self.nml['drv_in']['seq_timemgr_inparm']['start_ymd'] = int(self._run_end_date.strftime(date_fmt['cesm']))
            self._compute_run_dates()
            # - ML - Setting ydirini might not be needed, try without at some point
            self.nml['INPUT_IO']['gribin']['ydirini'] = self.nml['INPUT_IO']['ioctl']['ydir_restart_out']
            for gribout in self._get_gribouts():
                gribout['lwrite_const'] = False
            if not self.cosmo_only:
                self.nml['drv_in']['seq_infodata_inparm']['start_type'] = 'continue'
            self.write_open_nml()
            self._update_controller()

        return continue_run


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

def available(cls):
    if cls._target_machine is None:
        raise NotImplementedError("_target_machine class variable not set for Class {:s}".format(cls.__name__))
    else:
        available_cases[cls._target_machine] = cls
        return cls

@available
class daint_case(cc2_case):
    """Class defining a COSMO-CLM2 case on Piz Daint"""

    _target_machine='daint'
    _n_tasks_per_node = 12
    _default_install_dir = os.path.normpath(os.environ['SCRATCH'])


    def __init__(self, wall_time='24:00:00', account=None, partition=None,
                 shebang='#!/bin/bash', modules_opt='switch', pgi_version=None,
                 **base_case_args):

        self.wall_time = wall_time
        self.account = account
        self.modules_opt = modules_opt
        self.pgi_version = pgi_version
        self.shebang = shebang
        self.partition = partition
        cc2_case.__init__(self, **base_case_args)

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


    def _build_controller(self):

        logfile = '{:s}_{:s}-{:s}.out'.format(self.name,
                                              self._run_start_date.strftime(date_fmt['cesm']),
                                              self._run_end_date.strftime(date_fmt['cesm']))
        with open(os.path.join(self.path, 'controller'), mode='w') as script:
            # shebang
            script.write('{:s}\n\n'.format(self.shebang))

            # slurm options
            script.write('#SBATCH --constraint=gpu\n')
            script.write('#SBATCH --job-name={:s}\n'.format(self.name))
            script.write('#SBATCH --nodes={:d}\n'.format(self._n_nodes))
            script.write('#SBATCH --output={:s}\n'.format(logfile))
            script.write('#SBATCH --error={:s}\n'.format(logfile))
            script.write('#SBATCH --account={:s}\n'.format(self.account))
            script.write('#SBATCH --time={:s}\n'.format(self.wall_time))
            script.write('#SBATCH --gres=gpu:1\n')
            if self.partition is not None:
                script.write('#SBATCH --partition={:s}\n'.format(self.partition))
            script.write('\n')

            # environment variables
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
                script.write('export MPICH_G2G_PIPELINE=256\n')
                if self.cosmo_only:
                    script.write('export MPICH_RDMA_ENABLED_CUDA=1\n')
            script.write('\n')

            # Modules
            if self.modules_opt != 'none':
                # pgi programing environment
                if self.modules_opt == 'purge':
                    script.write('module purge\n')
                    script.write('module load PrgEnv-pgi\n')
                elif self.modules_opt == 'switch':
                    script.write('module switch PrgEnv-cray PrgEnv-pgi\n')
                # pgi version
                if self.pgi_version is not None:
                    script.write('module unload pgi\n')
                    script.write('module load pgi/{:s}\n'.format(self.pgi_version))

                # other modules
                script.write('module load daint-gpu\n')
                script.write('module load cray-netcdf\n')
                if self.gpu_mode:
                    script.write('module load craype-accel-nvidia60\n')    
                script.write('\n')

            # launch case
            script.write('cc2_control_case ./config.xml\n')


    def _update_controller(self):
        logfile = '{:s}_{:s}-{:s}.out'.format(self.name,
                                              self._run_start_date.strftime(date_fmt['cesm']),
                                              self._run_end_date.strftime(date_fmt['cesm']))
        rules = {'#SBATCH +--output=.*$': '#SBATCH --output={:s}'.format(logfile),
                 '#SBATCH +--error=.*$': '#SBATCH --error={:s}'.format(logfile)}
        with open(os.path.join(self.path, 'controller'), mode='r+') as f:
            content = f.read()
            for pattern, repl in rules.items():
                content = re.sub(pattern, repl, content, flags=re.MULTILINE)
            f.seek(0)
            f.write(content)
            f.truncate()


    def _update_xml_config(self, machine_node):
        ET.SubElement(machine_node, 'account').text = self.account
        ET.SubElement(machine_node, 'wall_time').text = self.wall_time
        ET.SubElement(machine_node, 'partition').text = self.partition
        ET.SubElement(machine_node, 'modules_opt').text = self.modules_opt
        ET.SubElement(machine_node, 'pgi_version').text = self.pgi_version
        ET.SubElement(machine_node, 'shebang').text = str(self.shebang)


    def _submit_func(self):
        check_call(['sbatch', 'controller'])


    def _run_func(self):
        # Determine run command
        if self.cosmo_only:
            if self.gpu_mode:
                run_cmd = 'srun -u --ntasks-per-node=1 -n {:d} {:s}'.format(self._n_nodes, self.COSMO_exe)
            else:
                run_cmd = 'srun -u -n {:d} {:s}'.format(self._n_nodes * self._n_tasks_per_node, self.COSMO_exe)
        else:
            self._build_proc_config()
            run_cmd = 'srun -u --multi-prog ./proc_config'

        # Run
        check_call(['module list'], shell=True)
        print("running " + run_cmd)
        sys.stdout.flush()
        check_call(run_cmd, shell=True)


    def _build_proc_config(self):

        # Build executable bash files
        f_path = os.path.join(self.path, 'cosmo.bash')
        with open(f_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("export MPICH_RDMA_ENABLED_CUDA={:1d}\n".format(self.gpu_mode))
            f.write("./{:s}".format(self.cos_exe))
        os.chmod(f_path, 0o755)
        f_path = os.path.join(self.path, 'cesm.bash')
        with open(f_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("export MPICH_RDMA_ENABLED_CUDA=0\n")
            f.write("./{:s}".format(self.cesm_exe))
        os.chmod(f_path, 0o755)

        # Build proc_config
        with open(os.path.join(self.path, 'proc_config'), mode='w') as f:
            if self.gpu_mode:
                N = self._n_tasks_per_node
                tasks = ",".join([str(k*N) for k in range(self._n_nodes)])
                f.write("{:s} ./cosmo.bash\n".format(tasks))
                tasks = ",".join(["{:d}-{:d}".format(k*N+1,(k+1)*N-1) for k in range(self._n_nodes)])
                f.write("{:s} ./cesm.bash".format(tasks))
            else:
                f.write('{:d}-{:d} ./cosmo.bash\n'.format(0, self._ncos-1))
                f.write('{:d}-{:d} ./cesm.bash'.format(self._ncos, self._ncos+self._ncesm-1))

@available
class mistral_case(cc2_case):
    """Class defining a COSMO-CLM2 case on Mistral"""

    _target_machine='mistral'
    _n_tasks_per_node = 24

    def __init__(self, wall_time='08:00:00', account=None, partition=None,
                 **base_case_args):
        self.wall_time = wall_time
        self.account = account
        self.partition = partition
        cc2_case.__init__(self, **base_case_args)
        if self.gpu_mode:
            raise NotImplementedError("gpu mode not implemented for " + self.__class__.__name__)


    def _build_proc_config(self):
        with open(os.path.join(self.path, 'proc_config'), mode='w') as f:
            f.write('{:d}-{:d} ./{:s}\n'.format(0, self._ncos-1, self.COSMO_exe))
            if not self.cosmo_only:
                f.write('{:d}-{:d} ./{:s}\n'.format(self._ncos, self._ncos+self._ncesm-1, self.CESM_exe))


    def _build_controller(self):
        logfile = '{:s}_{:s}-{:s}.out'.format(self.name,
                                              self._run_start_date.strftime(date_fmt['cesm']),
                                              self._run_end_date.strftime(date_fmt['cesm']))
        with open(os.path.join(self.path, 'controller'), mode='w') as script:
            # shebang
            script.write('#!/usr/bin/env bash\n')

            # slurm options
            script.write('#SBATCH --job-name={:s}\n'.format(self.name))
            script.write('#SBATCH --nodes={:d}\n'.format(self._n_nodes))
            script.write('#SBATCH --output={:s}\n'.format(logfile))
            script.write('#SBATCH --error={:s}\n'.format(logfile))
            script.write('#SBATCH --account={:s}\n'.format(self.account))
            script.write('#SBATCH --time={:s}\n'.format(self.wall_time))
            if self.partition is not None:
                script.write('#SBATCH --partition={:s}\n'.format(self.partition))
            script.write('\n')

            # environment variables
            script.write('export LD_LIBRARY_PATH=/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-parallel-openmpi2-intel14/lib/:/sw/rhel6-x64/netcdf/parallel_netcdf-1.6.1-openmpi2-intel14/lib\n')
            script.write('\n')
            script.write('# Set this to avoid segmentation faults\n')
            script.write('ulimit -s unlimited\n')
            script.write('ulimit -a\n')
            script.write('\n')
            script.write('export OMP_NUM_THREADS=1\n')
            script.write('\n')

            # launch case
            script.write('cc2_control_case ./config.xml\n')


    def _update_controller(self):
        logfile = '{:s}_{:s}-{:s}.out'.format(self.name,
                                              self._run_start_date.strftime(date_fmt['cesm']),
                                              self._run_end_date.strftime(date_fmt['cesm']))
        rules = {'#SBATCH +--output=.*$': '#SBATCH --output={:s}'.format(logfile),
                 '#SBATCH +--error=.*$': '#SBATCH --error={:s}'.format(logfile)}
        with open(os.path.join(self.path, 'controller'), mode='r+') as f:
            content = f.read()
            for pattern, repl in rules.items():
                content = re.sub(pattern, repl, content, flags=re.MULTILINE)
            f.seek(0)
            f.write(content)
            f.truncate()


    def _update_xml_config(self, config):
        ET.SubElement(config, 'machine').text = 'mistral'
        ET.SubElement(config, 'account').text = self.account
        ET.SubElement(config, 'wall_time').text = self.wall_time
        ET.SubElement(config, 'partition').text = self.partition


    def _submit_func(self):
        check_call(['sbatch', 'controller', './config.xml'])


    def _run_func(self):
        if self.cosmo_only:
            run_cmd = 'srun -u -n {:d} {:s}'.format(self._n_nodes * self._n_tasks_per_node, self.COSMO_exe)
        else:
            self._build_proc_config()
            run_cmd = 'srun -u --multi-prog ./proc_config'
        print("running " + run_cmd)
        sys.stdout.flush()
        check_call(run_cmd, shell=True)

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
        for name in self:
            self.write(name)

from __future__ import print_function
from .tools import date_fmt, add_time_from_str, COSMO_input_file_name, indent_xml
from subprocess import check_call, check_output
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

def available(cls):
    if cls._target_machine is None:
        raise NotImplementedError("_target_machine class variable not set for Class {:s}".format(cls.__name__))
    else:
        available_cases[cls._target_machine] = cls
        return cls

class cc2_case(object):
    """Base class defining a COSMO-CLM2 case"""

    _target_machine = None
    _n_tasks_per_node = None
    _default_install_dir = None
    _control_job = 'cc2_control_job'
    _run_job = 'cc2_run_job'
    _transfer_job = 'cc2_transfer_job'
    _archive_job = 'cc2_archive_job'
    _xml_config = 'cc2_config.xml'
    _transfer_list = 'cc2_transfer_list'
    NotImplementedMessage = "required method {:s} not implemented by class {:s}.\n" \
                            "Implement with a single pass statement if irrelevant to this machine."


    def __init__(self, name='COSMO_CLM2', install_dir=None, install=False,
                 cos_nml='./COSMO_nml', cos_in='./COSMO_input', cos_exe='./cosmo', cos_rst=None,
                 cesm_nml='./CESM_nml', cesm_in='./CESM_input', cesm_exe='./cesm.exe', cesm_rst=None,
                 oas_in='./OASIS_input', oas_nml='./OASIS_nml', archive_dir=None,
                 start_date=None, end_date=None, run_length=None,
                 ncosx=None, ncosy=None, ncosio=None, ncesm=None,
                 gpu_mode=False, dummy_day=True, cosmo_only=False, gen_oasis=False,
                 input_type='file', transfer_all=True,
                 archive_per_month=False, archive_compression='none', archive_cesm=True, archive_rm=False,
                 start_mode='startup', restart_date=None):

        # Basic init (no particular work required)
        self.name = name
        self.run_length = run_length
        self.gpu_mode = gpu_mode
        self.dummy_day = dummy_day
        self.cosmo_only = cosmo_only
        self.gen_oasis = gen_oasis
        self.cos_in = os.path.abspath(cos_in)
        self.install = install
        self.input_type = input_type
        self.transfer_all = transfer_all
        self.archive_dir = None if archive_dir is None else os.path.abspath(archive_dir)
        self.archive_per_month = archive_per_month
        self.archive_compression = archive_compression
        self.archive_cesm = archive_cesm
        self.archive_rm = archive_rm
        self.transfer_by_chunck = not self.transfer_all and self.input_type == 'file'
        self.start_mode = start_mode
        # Create namelists dictionnary
        self.nml = nmldict(self)
        # Set install_dir and path
        self.install_dir = install_dir
        # Install: transfer namelists, executables and input files
        if self.install:
            log = 'Setting up case {:s} in {:s}'.format(self.name, self._path)
            print(log + '\n' + '-' * len(log))
            self.install_case(cos_nml, cos_in, cos_exe, cos_rst, cesm_nml, cesm_in, cesm_exe, cesm_rst, oas_nml, oas_in)
            self.create_missing_dirs()
        self.cos_exe = cos_exe
        if not self.cosmo_only:
            self.cesm_exe = cesm_exe
        # Settings involving namelist changes
        self.start_date = start_date
        self.restart_date = restart_date
        self.end_date = end_date
        # - ML - Some of the following is useless for transfer and archive actions
        self._set_nml_start_parameters()
        self._compute_run_dates()   # defines _run_start_date, _run_end_date and _runtime (_end_date if needed)
        self._apply_run_dates()   # put dates and runtime in namelists objects (writing to file at the end)
        self._check_INPUT_IO()
        self._organize_tasks(ncosx, ncosy, ncosio, ncesm)
        # Finish install
        if self.install:
            self._build_run_job()
            if self.transfer_by_chunck and self._run_end_date < self.end_date:
                self._build_transfer_job()
            if self.archive_dir is not None:
                self._build_archive_job()
            self.to_xml()
            self.install_input()
        # Write modified namelists to file
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
    def start_date(self):
        return self._start_date
    @start_date.setter
    def start_date(self, start_date):
        if start_date is not None:
            self._start_date = datetime.strptime(start_date, date_fmt['in'])
            self.nml['INPUT_ORG']['runctl']['ydate_ini'] = self._start_date.strftime(date_fmt['cosmo'])
        elif 'ydate_ini' in self.nml['INPUT_ORG']['runctl']:
            self._start_date = datetime.strptime(self.nml['INPUT_ORG']['runctl']['ydate_ini'],
                                                 date_fmt['cosmo'])
        else:
            raise ValueError("ydate_ini has to be given in INPUT_ORG/runctl if no start_date is provided")

    @property
    def restart_date(self):
        return self._restart_date
    @restart_date.setter
    def restart_date(self, restart_date):
        if restart_date is not None:
            self._restart_date = datetime.strptime(restart_date, date_fmt['in'])
        else:
            self._restart_date = None

    @property
    def end_date(self):
        return self._end_date
    @end_date.setter
    def end_date(self, end_date):
        if end_date is not None:
            self._end_date = datetime.strptime(end_date, date_fmt['in'])
            self.nml['INPUT_ORG']['runctl']['ydate_end'] = self._end_date.strftime(date_fmt['cosmo'])
        elif 'ydate_end' in self.nml['INPUT_ORG']['runctl']:
            self._end_date = datetime.strptime(self.nml['INPUT_ORG']['runctl']['ydate_end'], date_fmt['cosmo'])
        else:
            self._end_date = None

    @property
    def cos_in_file_size(self):
        status_node = ET.parse(os.path.join(self.path, self._xml_config)).getroot().find('status')
        return int(status_node.find('cos_in_file_size').text)
    @cos_in_file_size.setter
    def cos_in_file_size(self, sz):
        tree = ET.parse(os.path.join(self.path, self._xml_config))
        tree.find('status').find('cos_in_file_size').text = str(sz)
        tree.write(os.path.join(self.path, self._xml_config), xml_declaration=True)

    @property
    def run_status(self):
        status_node = ET.parse(os.path.join(self.path, self._xml_config)).getroot().find('status')
        return status_node.find('run_status').text
    @run_status.setter
    def run_status(self, status):
        tree = ET.parse(os.path.join(self.path, self._xml_config))
        tree.find('status').find('run_status').text = status
        tree.write(os.path.join(self.path, self._xml_config), xml_declaration=True)

    @property
    def transfer_status(self):
        status_node = ET.parse(os.path.join(self.path, self._xml_config)).getroot().find('status')
        return status_node.find('transfer_status').text
    @transfer_status.setter
    def transfer_status(self, status):
        tree = ET.parse(os.path.join(self.path, self._xml_config))
        tree.find('status').find('transfer_status').text = status
        tree.write(os.path.join(self.path, self._xml_config), xml_declaration=True)


    def install_case(self, cos_nml, cos_in, cos_exe, cos_rst, cesm_nml, cesm_in, cesm_exe, cesm_rst, oas_nml, oas_in):

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

        if self.start_mode != 'startup':
            if 'ydir_restart_in' in self.nml['INPUT_IO']['ioctl']:
                cos_rst_nml = self.nml['INPUT_IO']['ioctl']['ydir_restart_in']
            elif 'ydir_restart' in  self.nml['INPUT_IO']['ioctl']:
                cos_rst_nml = self.nml['INPUT_IO']['ioctl']['ydir_restart']
            else:
                cos_rst_nml = ''
            self._mk_miss_path(cos_rst_nml)
            cos_rst_target = os.path.normpath(os.path.join(self.path, cos_rst_nml, os.path.basename(cos_rst)))
            check_call(['rsync', '-avL', cos_rst, cos_rst_target])

            if cos_rst.endswith('.gz'):
                check_call(['gunzip', cos_rst_target])
            elif cos_rst.endswith('.bz2'):
                check_call(['bzip2', cos_rst_target])

            if cesm_rst.endswith('.tar'):
                check_call(['tar', '-xvf', cesm_rst, '-C', self.path])
            elif cesm_rst.endswith(('.tgz', '.tar.gz')):
                check_call(['tar', '-zxvf', cesm_rst, '-C', self.path])
            elif cesm_rst.endswith(('.tbz', '.tar.bz2')):
                check_call(['tar', '-jxvf', cesm_rst, '-C', self.path])
            else:
                check_call(['rsync', '-avL', cesm_rst+'/', self.path])

        # Set case name in namelist
        self.nml['drv_in']['seq_infodata_inparm']['case_name'] = self.name


    def _set_nml_start_parameters(self):

        if self.start_mode == 'startup':
            self.nml['INPUT_ORG']['runctl']['hstart'] = 0
            if not self.cosmo_only:
                self.nml['drv_in']['seq_timemgr_inparm']['start_ymd'] = int(self.start_date.strftime(date_fmt['cesm']))
                self.nml['drv_in']['seq_infodata_inparm']['start_type'] = 'startup'
            if 'nrevsn' in self.nml['lnd_in']['clm_inparm']:
                del(self.nml['lnd_in']['clm_inparm']['nrevsn'])
        else:
            self.nml['INPUT_ORG']['runctl']['hstart'] = (self.restart_date - self.start_date).total_seconds() // 3600.0
            if not self.cosmo_only:
                self.nml['drv_in']['seq_timemgr_inparm']['start_ymd'] = int(self.restart_date.strftime(date_fmt['cesm']))
                if self.start_mode == 'continue':
                    self.nml['drv_in']['seq_infodata_inparm']['start_type'] = 'continue'
                    if 'nrevsn' in self.nml['lnd_in']['clm_inparm']:
                        del(self.nml['lnd_in']['clm_inparm']['nrevsn'])
                elif self.start_mode == 'restart':
                    self.nml['drv_in']['seq_infodata_inparm']['start_type'] = 'branch'
                    with open(os.path.join(self.path, 'rpointer.lnd')) as rpt:
                        self.nml['lnd_in']['clm_inparm']['nrevsn'] = rpt.readline().strip()


    def _cos_input_delta_ext(self):

        # Set time interval between 2 intput files
        delta = timedelta(hours=self.nml['INPUT_IO']['gribin']['hincbound'])
        # Set file extension
        ext = ''
        if 'yform_read' in self.nml['INPUT_IO']['ioctl']:
            if self.nml['INPUT_IO']['ioctl']['yform_read'] == 'ncdf':
                ext = '.nc'
        return delta, ext


    def build_transfer_list(self, start_date, end_date, initial=False):

        delta, ext = self._cos_input_delta_ext()

        # function to check and add file to transfer list or directly symlink
        def _check_add_file(root, date, file_list):
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
        with open(self._transfer_list, mode ='w') as t_list:
            if initial:
                _check_add_file('laf', start_date, t_list)
            cur_date = start_date
            while cur_date <= end_date:
                _check_add_file('lbfd', cur_date, t_list)
                cur_date += delta


    def transfer_input(self):

        if self.input_type == 'file':
            check_call(['rsync', '-avrL', '--files-from', self._transfer_list,
                        self.cos_in+'/', os.path.join(self.path,'COSMO_input')])


    def install_input(self):

        # Get cosmo lbf input file size
        _, ext = self._cos_input_delta_ext()
        file_name = COSMO_input_file_name('lbfd', self.start_date, ext)
        file_path = os.path.join(self.cos_in, file_name)
        self.cos_in_file_size = os.stat(file_path).st_size

        # Create COSMO_input directory if missing (essentially for symlinks)
        self._mk_miss_path('COSMO_input')

        # Transfer first chunck input files or all
        initial = self.start_mode == 'startup'
        if self.transfer_by_chunck and self._run_end_date < self.end_date:
            self.build_transfer_list(self._run_start_date, self._run_end_date, initial=initial)
        else:
            end_date = self.end_date + timedelta(days=1) if self.dummy_day else self.end_date
            self.build_transfer_list(self._run_start_date, end_date, initial=initial)
        self.transfer_input()
        os.remove(self._transfer_list)

        # Set transfer status
        self.transfer_status = 'complete'


    def _check_COSMO_input(self, start_date, end_date):

        delta, ext = self._cos_input_delta_ext()
        cur_date = start_date
        cos_in_file_size = self.cos_in_file_size # get from xml once for all
        while cur_date <= end_date:
            file_name = COSMO_input_file_name('lbfd', cur_date, ext)
            file_path = os.path.join(self.path, 'COSMO_input', file_name)
            if not os.path.exists(file_path):
                raise ValueError("COSMO input file {:s} missing".format(file_name))
            fs = os.stat(file_path).st_size
            if fs != cos_in_file_size:
                err_mess = "COSMO input file {:s} has byte size {:d} instead of {:d}"
                raise ValueError(err_mess.format(file_name, fs, cos_in_file_size))
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
        date_cosmo = self._start_date + timedelta(hours=INPUT_ORG['runctl']['hstart'])
        if not self.cosmo_only:
            date_cesm = datetime.strptime(str(drv_in['seq_timemgr_inparm']['start_ymd']), date_fmt['cesm'])
            if date_cosmo != date_cesm:
                raise ValueError("start dates are not identical in COSMO and CESM namelists")
        self._run_start_date = date_cosmo

        # Compute _runtime and _run_end_date (possibly _end_date)
        # -------------------------------------------------------
        if self._end_date is not None:
            if self._run_start_date >= self._end_date:
                raise ValueError("run sart date >= case end date")
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

        # Add a dummy day to the last chunk if required
        # ---------------------------------------------
        if self._run_end_date == self._end_date and self.dummy_day:
            one_day = timedelta(days=1)
            self._run_end_date += one_day
            self._runtime += one_day


    def _apply_run_dates(self):

        # Access to namelists
        INPUT_ORG = self.nml['INPUT_ORG']
        INPUT_IO = self.nml['INPUT_IO']
        if not self.cosmo_only:
            drv_in = self.nml['drv_in']

        # Compute times
        start_seconds = (self._run_start_date - self.start_date).total_seconds()
        dt = INPUT_ORG['runctl']['dt']
        runtime_seconds = self._runtime.total_seconds()
        runtime_hours = runtime_seconds // 3600

        hstart = start_seconds // 3600
        nstart = start_seconds // dt - 1
        hstop = hstart + runtime_hours
        nstop = (start_seconds + runtime_seconds) // dt


        # adapt INPUT_ORG
        if 'hstop' in INPUT_ORG['runctl']:
            del INPUT_ORG['runctl']['hstop']
        INPUT_ORG['runctl']['nstop'] = nstop - 1
        if 'hstop' in INPUT_ORG['runctl']:
            del INPUT_ORG['runctl']['hstop']

        # adapt INPUT_IO
        for gribout in self._get_gribouts():
            if 'hcomb' in gribout:
                gribout['hcomb'][0:2] = hstart, hstop
            elif 'ncomb' in gribout:
                gribout['ncomb'][0:2] = nstart, nstop
        if self._run_end_date > self.end_date:
            # ensure restart for the real end date, not after the dummy day
            INPUT_IO['ioctl']['nhour_restart'] = [int(hstop)-24, int(hstop)-24, 24]
        else:
            INPUT_IO['ioctl']['nhour_restart'] = [int(hstop), int(hstop), 24]

        if not self.cosmo_only:
            # adapt drv_in
            drv_in['seq_timemgr_inparm']['stop_n'] = int(runtime_seconds)
            if self._run_end_date > self.end_date:
                # ensure restart for the real end date, not after the dummy day
                drv_in['seq_timemgr_inparm']['restart_n'] = int(runtime_seconds)-86400
            else:
                drv_in['seq_timemgr_inparm']['restart_n'] = int(runtime_seconds)

            # adapt namcouple
            with open(os.path.join(self.path, 'namcouple_tmpl'), mode='r') as f:
                content = f.read()
            content = re.sub('_runtime_', str(int(runtime_seconds)), content)
            with open(os.path.join(self.path, 'namcouple'), mode='w') as f:
                f.write(content)


    def get_next_run_end_date(self):

        next_end_date = min(add_time_from_str(self._run_end_date, self.run_length),
                            self.end_date)
        if next_end_date == self.end_date and self.dummy_day:
            next_end_date += timedelta(days=1)

        return next_end_date


    def _check_INPUT_IO(self):

        # Make sure COSMO input and initial files are looked for in the COSMO_input folder
        self.nml['INPUT_IO']['gribin']['ydirini'] = 'COSMO_input'
        self.nml['INPUT_IO']['gribin']['ydirbd'] = 'COSMO_input'

        # Only keep gribout blocks that fit within runtime
        # (essentially to avoid crash for short tests)
        runtime_seconds = self._runtime.total_seconds()
        runtime_hours = runtime_seconds // 3600.0
        runtime_nstep = runtime_seconds // self.nml['INPUT_ORG']['runctl']['dt']
        gribouts_out = []
        gribouts_in = self._get_gribouts()
        for gribout in gribouts_in:
            if 'hcomb' in gribout:
                if runtime_hours >= gribout['hcomb'][2]:
                    gribouts_out.append(gribout)
            elif 'ncomb' in gribout:
                if runtime_nstep >= gribout['ncomb'][2]:
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


    def create_missing_dirs(self):

        # COSMO
        # -----
        # input
        self._mk_miss_path(self.nml['INPUT_IO']['gribin']['ydirini'])
        self._mk_miss_path(self.nml['INPUT_IO']['gribin']['ydirbd'])
        # output
        for gribout in self._get_gribouts():
            self._mk_miss_path(gribout['ydir'])
        if 'ydir_restart' in self.nml['INPUT_IO']['ioctl']:
            self._mk_miss_path(self.nml['INPUT_IO']['ioctl']['ydir_restart'])
        if 'ydir_restart_in' in self.nml['INPUT_IO']['ioctl']:
            self._mk_miss_path(self.nml['INPUT_IO']['ioctl']['ydir_restart_in'])
        if 'ydir_restart_out' in self.nml['INPUT_IO']['ioctl']:
            self._mk_miss_path(self.nml['INPUT_IO']['ioctl']['ydir_restart_out'])

        # CESM
        # ----
        if not self.cosmo_only:
            # timing
            # remove if exists before creating
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


    def to_xml(self):

        config_node = ET.Element('config')
        tree = ET.ElementTree(config_node)
        ET.SubElement(config_node, 'machine').text = self._target_machine

        main_node = ET.SubElement(config_node, 'main')
        ET.SubElement(main_node, 'name').text = self.name
        ET.SubElement(main_node, 'install_dir').text = self.install_dir
        ET.SubElement(main_node, 'cosmo_only', type='py_eval').text = str(self.cosmo_only)
        ET.SubElement(main_node, 'gen_oasis', type='py_eval').text = str(self.gen_oasis)
        ET.SubElement(main_node, 'run_length').text = self.run_length
        ET.SubElement(main_node, 'cos_exe').text = self.cos_exe
        if not self.cosmo_only:
            ET.SubElement(main_node, 'cesm_exe').text = self.cesm_exe
        ET.SubElement(main_node, 'cos_in').text = self.cos_in
        ET.SubElement(main_node, 'archive_dir').text = self.archive_dir
        ET.SubElement(main_node, 'gpu_mode', type='py_eval').text = str(self.gpu_mode)
        ET.SubElement(main_node, 'dummy_day', type='py_eval').text = str(self.dummy_day)
        ET.SubElement(main_node, 'transfer_all', type='py_eval').text = str(self.transfer_all)
        ET.SubElement(main_node, 'start_mode').text = self.start_mode
        node = ET.SubElement(main_node, 'restart_date')
        if self.start_mode != 'startup':
            node.text = self.restart_date.strftime(date_fmt['in'])

        status_node = ET.SubElement(config_node, 'status')
        ET.SubElement(status_node, 'run_status')
        ET.SubElement(status_node, 'transfer_status')
        ET.SubElement(status_node, 'cos_in_file_size')

        # - ML - Could be usefull in case machine specific arguments need to be stored one day.
        #        This isn't the case as of now
        ET.SubElement(config_node, self._target_machine)

        indent_xml(config_node)

        tree.write(os.path.join(self.path, self._xml_config), xml_declaration=True)


    def set_next_run(self):

        if self._run_end_date < self._end_date:
            tree = ET.parse(os.path.join(self.path, self._xml_config))
            main_node = tree.find('main')
            main_node.find('start_mode').text = 'continue'
            main_node.find('restart_date').text = self._run_end_date.strftime(date_fmt['in'])
            tree.write(os.path.join(self.path, self._xml_config), xml_declaration=True)


    def submit_run(self):

        cwd = os.getcwd()
        os.chdir(self.path)

        self._submit_run_cmd(self._run_start_date, self._run_end_date)

        os.chdir(cwd)


    def submit_next_run(self):

        cwd = os.getcwd()
        os.chdir(self.path)

        self._submit_run_cmd(self._run_end_date, self.get_next_run_end_date())

        os.chdir(cwd)


    def submit_next_transfer(self):

        cwd = os.getcwd()
        os.chdir(self.path)

        next_end_date = self.get_next_run_end_date()
        self.build_transfer_list(self._run_end_date, next_end_date)
        self._submit_transfer_cmd(self._run_end_date, next_end_date)

        os.chdir(cwd)


    def submit_archive(self):

        cwd = os.getcwd()
        os.chdir(self.path)

        self._submit_archive_cmd()

        os.chdir(cwd)


    def run(self):

        # Monitor time
        start_time = time.time()

        # Clean workdir
        file_list = glob('YU*') + glob('debug*') + glob('core*') + glob('nout.*') + glob('*.timers_*')
        for f in file_list:
            os.remove(f)

        # Check presence and size of input files for current chunk
        self._check_COSMO_input(self._run_start_date, self._run_end_date)

        # Run
        self._run_fun()

        # Monitor time
        elapsed = time.time() - start_time
        print("\nCase {name:s} ran in {elapsed:.2f}\n".format(name=self.name, elapsed=elapsed))


    def _build_run_job(self):
        """Place holder for _build_run_job method to be implemented by machine specific classes."""

        raise NotImplementedError(self.NotImplementedMessage.format('_build_run_job(self)', self.__class__.__name__))


    def _build_transfer_job(self):
        """Place holder for _build_transfer_job method to be implemented by machine specific classes."""

        raise NotImplementedError(self.NotImplementedMessage.format('_build_transfer_job(self)', self.__class__.__name__))


    def _build_archive_job(self):
        """Place holder for _build_archive_job method to be implemented by machine specific classes."""

        raise NotImplementedError(self.NotImplementedMessage.format('_build_archive_job(self)', self.__class__.__name__))


    def _run_fun(self):
        """Place holder for _run_fun method to be implemented by machine specific classes."""

        raise NotImplementedError(self.NotImplementedMessage.format('_run_fun(self)', self.__class__.__name__))


    def _submit_run_cmd(self):
        """Place holder for _submit_run_cmd method to be implemented by machine specific classes."""

        raise NotImplementedError(self.NotImplementedMessage.format('_submit_run_cmd(self)', self.__class__.__name__))


    def _submit_transfer_cmd(self):
        """Place holder for _submit_transfer_cmd method to be implemented by machine specific classes."""

        raise NotImplementedError(self.NotImplementedMessage.format('_submit_transfer_cmd(self)', self.__class__.__name__))


    def _submit_archive_cmd(self):
        """Place holder for _submit_archive_cmd method to be implemented by machine specific classes."""

        raise NotImplementedError(self.NotImplementedMessage.format('_submit_archive_cmd(self)', self.__class__.__name__))

@available
class daint_case(cc2_case):
    """Class defining a COSMO-CLM2 case on Piz Daint"""

    _target_machine='daint'
    _n_tasks_per_node = 12
    _default_install_dir = os.path.normpath(os.environ['SCRATCH'])
    _archive_rst_job = 'cc2_archive_rst_job'


    def __init__(self, run_time='24:00:00', account=None, partition=None,
                 shebang='#!/bin/bash', modules_opt='switch', pgi_version=None,
                 transfer_time='02:00:00', archive_time='03:00:00', **base_case_args):

        self.run_time = run_time
        self.transfer_time = transfer_time
        self.archive_time = archive_time
        self.account = account
        self.modules_opt = modules_opt
        self.pgi_version = pgi_version
        self.shebang = shebang
        self.partition = partition
        cc2_case.__init__(self, **base_case_args)
        if self.install:
            if not self.cosmo_only:
                self._build_proc_config()
            self.update_xml_config()


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


    def update_xml_config(self):
        tree = ET.parse(os.path.join(self.path, self._xml_config))
        daint_node = tree.find('daint')
        ET.SubElement(daint_node, 'archive_per_month', type='py_eval').text = str(self.archive_per_month)
        indent_xml(tree.getroot())
        tree.write(os.path.join(self.path, self._xml_config), xml_declaration=True)


    def _build_run_job(self):

        # Header
        # ------
        script_str = '{:s}\n\n'.format(self.shebang)
        script_str += '#SBATCH --constraint=gpu\n'
        script_str += '#SBATCH --job-name={:s}\n'.format(self.name)
        script_str += '#SBATCH --nodes={:d}\n'.format(self._n_nodes)
        script_str += '#SBATCH --account={:s}\n'.format(self.account)
        script_str += '#SBATCH --time={:s}\n'.format(self.run_time)
        script_str += '#SBATCH --gres=gpu:1\n'
        if self.partition is not None:
            script_str += '#SBATCH --partition={:s}\n'.format(self.partition)

        # environment variables
        # ---------------------
        script_str +='''
export MALLOC_MMAP_MAX_=0
export MALLOC_TRIM_THRESHOLD_=536870912

# Set this to avoid segmentation faults
ulimit -s unlimited
ulimit -a

export OMP_NUM_THREADS=1
'''
        if self.gpu_mode:
            script_str += '''
# Use for gpu mode
export MV2_ENABLE_AFFINITY=0
export MV2_USE_CUDA=1
export MPICH_G2G_PIPELINE=256
'''
            if self.cosmo_only:
                script_str += 'export MPICH_RDMA_ENABLED_CUDA=1\n'

        # Modules
        # -------
        script_str += '\n'
        if self.modules_opt != 'none':
            # pgi programing environment
            if self.modules_opt == 'purge':
                script_str += 'module purge\n'
                script_str += 'module load PrgEnv-pgi\n'
            elif self.modules_opt == 'switch':
                script_str += 'module switch PrgEnv-cray PrgEnv-pgi\n'
            # pgi version
            if self.pgi_version is not None:
                script_str += 'module unload pgi\n'
                script_str += 'module load pgi/{:s}\n'.format(self.pgi_version)

            # other modules
            script_str += 'module load daint-gpu\n'
            script_str += 'module load cray-netcdf\n'
            if self.gpu_mode:
                script_str += 'module load craype-accel-nvidia60\n'
            script_str += '\n'

        # launch case
        # -----------
        script_str += 'cc2_control_case ./{:s}'.format(self._xml_config)

        # Write to file
        # -------------
        with open(os.path.join(self.path, self._run_job), mode='w') as script:
            script.write(script_str)


    def _submit_run_cmd(self, d1, d2):

        d1_str = d1.strftime(date_fmt['cesm'])
        d2_str = d2.strftime(date_fmt['cesm'])
        logfile = '{:s}_{:s}-{:s}.out'.format(self.name, d1_str, d2_str)

        cmd_tmpl = 'sbatch --output={log:s} --error={log:s} {job:s}'
        cmd = cmd_tmpl.format(log=logfile, job=self._run_job)
        print("submitting run with check_call('" + cmd + "', shell=True)")
        check_call(cmd, shell=True)


    def _build_transfer_job(self):

        # Header
        # ------
        script_str = '#!/bin/bash -l\n\n'
        script_str += '#SBATCH --partition=xfer\n'
        script_str += '#SBATCH --ntasks=1\n'
        script_str += '#SBATCH --time={:s}\n'.format(self.transfer_time)
        script_str += '#SBATCH --job-name=cc2_transfer\n\n'

        # Case dependent variables
        # ------------------------
        script_str += '# Case dependent variables\n'
        script_str += '# ------------------------\n'
        script_str += 'run_job={:s}\n'.format(self._run_job)
        script_str += 'xml_config={:s}\n'.format(self._xml_config)
        script_str += 'cos_in_origin={:s}\n'.format(self.cos_in)
        script_str += 'cos_in_target={:s}\n'.format(os.path.join(self.path,'COSMO_input'))
        script_str += 'transfer_list={:s}\n'.format(self._transfer_list)

        # Main script
        # -----------
        # Use sed commands to handle case status as python isn't available on the xfer queue.
        # Otherwise just use cc2_control --action=transfer
        script_str += '''
# Functions to get and set case status
# ------------------------------------
get_status(){
    sed -n 's@\s*<'$1'_status>\(.*\)</'$1'_status>@\\1@p' ${xml_config}
}
set_status(){
    sed -i 's@\(\s*<'$1'_status>\).*\(</'$1'_status>\)@\\1'$2'\\2@' ${xml_config}
}

# Set transfer status
# -------------------
set_status "transfer" "transferring"

# Transfer
# --------
rsync -avrL --files-from ${transfer_list} ${cos_in_origin} ${cos_in_target}

# Submit next run
# ---------------
if [[ $(get_status "run") == "complete" ]]; then
    set_status "run" "submitted"
    unset SLURM_MEM_PER_CPU
    sbatch --output $1 --error $1 ${run_job}
fi

# Set transfer status
# -------------------
set_status "transfer" "complete"'''

        # Write to file
        # -------------
        with open(os.path.join(self.path, self._transfer_job), mode='w') as script:
            script.write(script_str)


    def _submit_transfer_cmd(self, d1, d2):

        d1_str = d1.strftime(date_fmt['cesm'])
        d2_str = d2.strftime(date_fmt['cesm'])
        logfile = 'transfer_{:s}-{:s}.out'.format(d1_str, d2_str)
        run_log = '{:s}_{:s}-{:s}.out'.format(self.name, d1_str, d2_str)

        cmd_tmpl = 'sbatch --output={log:s} --error={log:s} {job:s} {run_log:s}'
        cmd = cmd_tmpl.format(log=logfile, job=self._transfer_job, run_log=run_log)
        print("submitting transfer with check_call('" + cmd + "', shell=True)")
        check_call(cmd, shell=True)


    def _build_archive_job(self):

        # Build archive job for output files
        # ==================================

        # Header
        # ------
        script_str = '#!/bin/bash -l\n\n'
        script_str += '#SBATCH --partition=xfer\n'
        script_str += '#SBATCH --ntasks=1\n'
        script_str += '#SBATCH --time={:s}\n'.format(self.archive_time)
        script_str += '#SBATCH --job-name=cc2_archive\n\n'

        # Case dependent variables
        # ------------------------
        script_str += '# Case dependent variables\n'
        script_str += '# ------------------------\n'
        script_str += 'CASE_NAME={:s}\n'.format(self.name)
        script_str += 'archive_dir={:s}\n'.format(os.path.join(self.archive_dir, self.name))
        script_str += 'archive_cesm={:s}\n'.format('true' if self.archive_cesm else 'false')
        # COSMO output streams
        stream_list = ['"{:s}"'.format(os.path.normpath(gribout['ydir'])) for gribout in self._get_gribouts()]
        script_str += 'COSMO_gribouts=({:s})\n'.format(' '.join(stream_list))
        # CESM output streams
        stream_list = ['"h0"']
        for k in range(2,7):
           if 'hist_fincl{:d}'.format(k) in self.nml['lnd_in']['clm_inparm']:
               stream_list += ['"h{:d}"'.format(k-1)]
        script_str += 'CESM_hh=({:s})\n'.format(' '.join(stream_list))
        # Archiving options
        script_str += 'remove_originals={:s}\n'.format('true' if self.archive_rm else 'false')
        script_str += 'compression={:s}\n'.format(self.archive_compression)

        # Main script
        # -----------
        script_str += '''
# Extract date components
# -----------------------
YS="${1:0:4}"
MS=$(echo "${1:4:2}" | sed 's/^0*//')
YE="${2:0:4}"
ME=$(echo "${2:4:2}" | sed 's/^0*//')

# Archiving options
# -----------------
if [[ ${compression} == "gzip" ]]; then
    tar_ext="tar.gz"
    tar_opt="zcf"
elif [[ ${compression} == "bzip2" ]]; then
    tar_ext="tar.bz2"
    tar_opt="jcf"
else
    tar_ext="tar"
    tar_opt="cf"
fi

if [[ ${remove_originals} == "true" ]]; then
    tar_rm='--remove-files'
else
    tar_rm=''
fi

# Proceed to archiving
# --------------------
mkdir -p ${archive_dir}/COSMO_output
mkdir -p ${archive_dir}/CESM_output

for ((YYYY=YS; YYYY<=YE; YYYY++)); do
    echo "treating year ${YYYY}"
    if [[ $YYYY == $YS ]]; then m1=$MS; else m1=01; fi
    if [[ $YYYY == $YE ]]; then m2=$ME; else m2=12; fi
    for ((m=m1; m<=m2; m++)); do
        echo "    treating month ${m}"
        # Archive COSMO output
        if ((${#COSMO_gribouts[@]} > 0)); then
            YYYYMM=${YYYY}$(printf "%02d" ${m})
            YYYYMMp1=$((YYYY + m/12))$(printf "%02d" $((m%12+1)))
            for gribout in ${COSMO_gribouts[@]}; do
                echo "        handling COSMO stream ${gribout}"
                gribname=$(basename ${gribout})
                arch_name=lffd${YYYYMM}.${tar_ext}
                files=$(find ${gribout} \( \( -name "lffd${YYYYMM}"'*' -and -not -name "lffd${YYYYMM}0100"'*' \)\\
                             -or -name "lffd${YYYYMMp1}0100"'*' \) -printf '%f\\n' | sort)
                if (( ${#files} > 0 )); then
                    echo "            preparing ${arch_name}"
                    tar -${tar_opt} ${arch_name} -C ${gribout} ${files} ${tar_rm}
                    mkdir -p ${archive_dir}/COSMO_output/${gribname}
                    echo "            sending ${arch_name} to archive directory"
                    rsync -ar ${arch_name} ${archive_dir}/COSMO_output/${gribname} --remove-source-files
                fi
            done
        fi
        # Archive CESM output
        if [[ ${#CESM_hh[@]} > 0 && ${archive_cesm} == "true" ]]; then
            YYYYMM=${YYYY}-$(printf "%02d" ${m})
            for hh in ${CESM_hh[@]}; do
                echo "        handling CESM stream ${hh}"
                arch_name=${CASE_NAME}.clm2.${hh}.${YYYYMM}.${tar_ext}
                files=$(find . -name "${CASE_NAME}.clm2.${hh}.${YYYYMM}"'*' -printf '%f\\n' | sort)
                if (( ${#files} > 0 )); then
                    echo "            preparing ${arch_name}"
                    tar -${tar_opt} ${arch_name} ${files} ${tar_rm}
                    mkdir -p ${archive_dir}/CESM_output/${hh}
                    echo "            sending ${arch_name} to archive directory"
                    rsync -ar ${arch_name} ${archive_dir}/CESM_output/${hh} --remove-source-files
                fi
            done
        fi
    done
done'''

        # Write to file
        # -------------
        with open(os.path.join(self.path, self._archive_job), mode='w') as script:
            script.write(script_str)

        # Build archive job for restart files
        # ===================================

        # Header
        # ------
        script_str = '#!/bin/bash -l\n\n'
        script_str += '#SBATCH --partition=xfer\n'
        script_str += '#SBATCH --ntasks=1\n'
        script_str += '#SBATCH --time={:s}\n'.format(self.archive_time)
        script_str += '#SBATCH --job-name=cc2_archive_rst\n\n'

        # Case dependent variables
        # ------------------------
        script_str += '# Case dependent variables\n'
        script_str += 'CASE_NAME={:s}\n'.format(self.name)
        script_str += 'archive_dir={:s}\n'.format(os.path.join(self.archive_dir, self.name))
        if 'ydir_restart' in self.nml['INPUT_IO']['ioctl']:
            COSMO_restart_dir = self.nml['INPUT_IO']['ioctl']['ydir_restart']
        elif 'ydir_restart_in' in self.nml['INPUT_IO']['ioctl']:
            COSMO_restart_dir = self.nml['INPUT_IO']['ioctl']['ydir_restart_in']
        else:
            COSMO_restart_dir = '.'
        script_str += 'COSMO_restart_dir={:s}\n\n'.format(COSMO_restart_dir)
        script_str += 'compression={:s}\n'.format(self.archive_compression)

        # Main script
        # -----------
        script_str += '''
# Extract date components
# -----------------------
YYYY="${1:0:4}"
MM="${1:4:2}"
DD="${1:6:2}"
HH="${1:8:2}"


# Archive COSMO restart file
# --------------------------
echo "handling COSMO restart file"
cd ${COSMO_restart_dir}
rst_file=lrfd${YYYY}${MM}${DD}${HH}o

# Compress rstart file
echo "    compressing ${rst_file} if needed"
if [[ ${compression} == "none" ]]; then
    transfer_file=${rst_file}
elif [[ ${compression} == "gzip" ]]; then
    transfer_file=${rst_file}.gz
    gzip < ${rst_file} > ${transfer_file}
elif [[ ${compression} == "bzip2" ]]; then
    transfer_file=${rst_file}.bz2
    bzip2 < ${rst_file} > ${transfer_file}
fi

# Transfer file
target_dir=${archive_dir}/COSMO_output/restart
mkdir -p ${target_dir}
echo "    sending ${transfer_file} to ${target_dir}"
rsync -ar ${transfer_file} ${target_dir}
if [[ $? == 0 && ${compression} != "none" ]]; then rm -f ${transfer_file}; fi

cd - > /dev/null


# Archive CESM restart files
# --------------------------
echo "handling CESM restart files"
date=${YYYY}-${MM}-${DD}-$(printf "%05d" $(( ${HH} * 3600 )) )
tar_name=restart_${date}.tar

# Generate rpointer files
echo "    generating rpointer files"
tmp_dir=tmp_archive_rst_${date}
mkdir -p ${tmp_dir}

printf "%s\n%s" ${CASE_NAME}.datm.r.${date}.nc ${CASE_NAME}.datm.rs1.${date}.bin > ${tmp_dir}/rpointer.atm
echo "${CASE_NAME}.clm2.r.${date}.nc" > ${tmp_dir}/rpointer.lnd
echo "${CASE_NAME}.cpl.r.${date}.nc" > ${tmp_dir}/rpointer.drv

# Start archive with rpointer files
echo "    building ${tar_name}"
tar -cf ${tar_name} -C ${tmp_dir} rpointer.*

rm -rf ${tmp_dir}

# Add CESM restart files
files=$( find \( -path "./${CASE_NAME}.clm2.rh[0-6].${date}.nc"\\
              -or -path "./${CASE_NAME}.clm2.r.${date}.nc"\\
              -or -path "./${CASE_NAME}.cpl.r.${date}.nc"\\
              -or -path "./${CASE_NAME}.datm.rs1.${date}.bin" \)\\
              -printf '%f\\n' | sort )
tar -rf ${tar_name} ${files}

# Compress archive
echo "    compressing ${tar_name} if needed"
if [[ ${compression} == "gzip" ]]; then
    transfer_name=${tar_name}.gz
    gzip < ${tar_name} > ${transfer_name} && rm ${tar_name}
elif [[ ${compression} == "bzip2" ]]; then
    transfer_name=${tar_name}.bz2
    bzip2 < ${tar_name} > ${transfer_name} && rm ${tar_name}
else
    transfer_name=${tar_name}
fi

# Transer archive
target_dir=${archive_dir}/CESM_output/restart
mkdir -p ${target_dir}
echo "    transferring ${transfer_name} to ${target_dir}"
rsync -ar ${transfer_name} ${target_dir} --remove-source-files'''

        # Write to file
        # -------------
        with open(os.path.join(self.path, self._archive_rst_job), mode='w') as script:
            script.write(script_str)


    def _submit_archive_cmd(self):

        # Submit output archive job(s)
        # ============================
        def _assemble_cmd_and_submit(d1, d2):
            d1_str, d2_str = d1.strftime('%Y%m'), d2.strftime('%Y%m')
            cmd_tmpl = 'sbatch --output={log:s} --error={log:s} {job:s} {d1:s} {d2:s}'
            logfile = '{:s}_{:s}-{:s}.out'.format('archive', d1_str, d2_str)
            cmd = cmd_tmpl.format(job=self._archive_job, d1=d1_str, d2=d2_str, log=logfile)
            print("submitting archive with check_call('" + cmd + "', shell=True)")
            check_call(cmd, shell=True)

        # Shift archiving period vs run period by 1 month except first and last chunks
        # (last COSMO output files written at the sart of next chunk)
        start_archive_date = max(self.start_date, add_time_from_str(self._run_start_date, '-1m'))
        if self._run_end_date >= self.end_date:
            end_archive_date = self.end_date
        else:
            end_archive_date = max(self.start_date, add_time_from_str(self._run_end_date, '-1m'))

        if self.archive_per_month:
            # Submit one archive job per month
            cur_start_date = start_archive_date
            while cur_start_date < end_archive_date:
                cur_month_beg = datetime(cur_start_date.year, cur_start_date.month, 1)
                tmp_date = add_time_from_str(cur_start_date, '1m')
                cur_month_end = datetime(tmp_date.year, tmp_date.month, 1)
                # Check that the we have to proceed to archiving for that period
                if cur_month_end <= end_archive_date or end_archive_date >= self.end_date:
                    # Determine date for archive job and submit
                    _assemble_cmd_and_submit(cur_month_beg, cur_month_beg)
                # Shift dates
                cur_start_date = cur_month_end
        else:
            # Submit the whole archiving period as one job
            proceed = False
            # Determine dates for archive job
            if end_archive_date >= self.end_date:
                d1 = datetime(start_archive_date.year, start_archive_date.month, 1)
                d2 = datetime(end_archive_date.year, end_archive_date.month, 1)
                if d2 == end_archive_date:
                    d2 = add_time_from_str(end_archive_date, '-1m')
                proceed = True
            else:
                tmp_date = add_time_from_str(start_archive_date, '1m')
                first_month_end = datetime(tmp_date.year, tmp_date.month, 1)
                if end_archive_date >= first_month_end:
                    d1 = datetime(start_archive_date.year, start_archive_date.month, 1)
                    tmp_date = add_time_from_str(end_archive_date, '-1m')
                    d2 = datetime(tmp_date.year, tmp_date.month, 1)
                    proceed = True
            if proceed:
                # Submit archive job
                _assemble_cmd_and_submit(d1, d2)

        # Submit restart archive job
        # ==========================
        end_date = min(self._run_end_date, self.end_date)
        end_date_str = end_date.strftime(date_fmt['cosmo'])
        logfile = '{:s}_{:s}.out'.format('archive_restart', end_date_str)
        cmd_tmpl = 'sbatch --output={log:s} --error={log:s} {job:s} {d:s}'
        cmd = cmd_tmpl.format(job=self._archive_rst_job, d=end_date_str, log=logfile)
        print("submitting restart archive with check_call('" + cmd + "', shell=True)")
        check_call(cmd, shell=True)


    def _run_fun(self):
        # Determine run command
        if self.cosmo_only:
            if self.gpu_mode:
                run_cmd = 'srun -u --ntasks-per-node=1 -n {:d} {:s}'.format(self._n_nodes, self.cos_exe)
            else:
                run_cmd = 'srun -u -n {:d} {:s}'.format(self._n_nodes * self._n_tasks_per_node, self.cos_exe)
        else:
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
            if self.gpu_mode: 
                f.write("export MPICH_RDMA_ENABLED_CUDA=1\n")
            f.write("./{:s}".format(self.cos_exe))
        os.chmod(f_path, 0o755)
        f_path = os.path.join(self.path, 'cesm.bash')
        with open(f_path, 'w') as f:
            f.write("#!/bin/bash\n")
            if self.gpu_mode:
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


    def __init__(self, run_time='08:00:00', account=None, partition=None,
                 transfer_time='02:00:00', archive_time='03:00:00', **base_case_args):

        self.run_time = run_time
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


    def _build_run_job(self):

        # shebang
        script_str = '#!/usr/bin/env bash\n\n'

        # slurm options
        scripte_str +='#SBATCH --job-name={:s}\n'.format(self.name)
        scripte_str +='#SBATCH --nodes={:d}\n'.format(self._n_nodes)
        scripte_str +='#SBATCH --account={:s}\n'.format(self.account)
        scripte_str +='#SBATCH --time={:s}\n'.format(self.run_time)
        if self.partition is not None:
            scripte_str += '#SBATCH --partition={:s}\n'.format(self.partition)

        # environment variables
        scripte_str += 'export LD_LIBRARY_PATH=/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-parallel-openmpi2-intel14/lib/:/sw/rhel6-x64/netcdf/parallel_netcdf-1.6.1-openmpi2-intel14/lib\n\n'
        script_str += '# Set this to avoid segmentation faults\n'
        script_str += 'ulimit -s unlimited\n'
        script_str += 'ulimit -a\n\n'
        script_str += 'export OMP_NUM_THREADS=1\n\n'

        # launch case
        scripte_str += 'cc2_control_case ./{:s}\n'.format(self._xml_config)

        with open(os.path.join(self.path, self._run_job), mode='w') as script:
            script.write(script_str)


    def _submit_run_cmd(self, d1, d2):

        d1_str = d1.strftime(date_fmt['cesm'])
        d2_str = d2.strftime(date_fmt['cesm'])
        logfile = '{:s}_{:s}-{:s}.out'.format(self.name, d1_str, d2_str)

        cmd_tmpl = 'sbatch --output={log:s} --error={log:s} {job:s}'
        cmd = cmd_tmpl.format(log=logfile, job=self._run_job)
        print("submitting run with check_call('" + cmd + "', shell=True)")
        check_call(cmd, shell=True)


    def _run_fun(self):
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

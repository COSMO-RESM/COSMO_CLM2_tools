from __future__ import print_function
from .base_case import base_case
from .date_formats import date_fmt_in, date_fmt_cosmo, date_fmt_cesm
from subprocess import check_call
import os
import re
import xml.etree.ElementTree as ET
import sys


class mistral_case(base_case):
    """Class defining a COSMO-CLM2 case on Mistral"""

    _n_tasks_per_node = 24

    def __init__(self, wall_time='08:00:00', account=None, partition=None,
                 **base_case_args):
        self.wall_time = wall_time
        self.account = account
        self.partition = partition
        base_case.__init__(self, **base_case_args)
        if self.gpu_mode:
            raise NotImplementedError("gpu mode not implemented for " + self.__class__.__name__)


    def _build_proc_config(self):
        with open(os.path.join(self.path, 'proc_config'), mode='w') as f:
            f.write('{:d}-{:d} ./{:s}\n'.format(0, self._ncos-1, self.COSMO_exe))
            if not self.cosmo_only:
                f.write('{:d}-{:d} ./{:s}\n'.format(self._ncos, self._ncos+self._ncesm-1, self.CESM_exe))


    def _build_controller(self):
        logfile = '{:s}_{:s}-{:s}.out'.format(self.name,
                                              self._run_start_date.strftime(date_fmt_cesm),
                                              self._run_end_date.strftime(date_fmt_cesm))
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

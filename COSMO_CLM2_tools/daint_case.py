from __future__ import print_function
from .base_case import base_case
from .date_formats import date_fmt_in, date_fmt_cosmo, date_fmt_cesm
from subprocess import check_call
import os
import re
import xml.etree.ElementTree as ET
import sys


class daint_case(base_case):
    """Class defining a COSMO-CLM2 case on Piz Daint"""

    _n_tasks_per_node = 12

    def __init__(self, wall_time='24:00:00', account=None, partition=None,
                 login_shell=True, modules_opt='switch', pgi_version=None,
                 **base_case_args):
        self.wall_time = wall_time
        self.account = account
        self.modules_opt = modules_opt
        self.pgi_version = pgi_version
        self.login_shell = login_shell
        self.partition = partition
        base_case.__init__(self, **base_case_args)


    def _build_proc_config(self):
        with open(os.path.join(self.path, 'proc_config'), mode='w') as f:
            if self.gpu_mode:
                N = self._n_tasks_per_node
                line = ",".join([str(k*N) for k in range(self._n_nodes)])
                f.write("{:s} ./{:s}\n".format(line, self.COSMO_exe))
                if not self.cosmo_only:
                    line = ",".join(["{:d}-{:d}".format(k*N+1,(k+1)*N-1) for k in range(self._n_nodes)])
                    f.write("{:s} ./{:s}\n".format(line, self.CESM_exe))
            else:
                f.write('{:d}-{:d} ./{:s}\n'.format(0, self._ncos-1, self.COSMO_exe))
                if not self.cosmo_only:
                    f.write('{:d}-{:d} ./{:s}\n'.format(self._ncos, self._ncos+self._ncesm-1, self.CESM_exe))


    def _build_controller(self):
        logfile = '{:s}_{:s}-{:s}.out'.format(self.name,
                                              self._run_start_date.strftime(date_fmt_cesm),
                                              self._run_end_date.strftime(date_fmt_cesm))
        with open(os.path.join(self.path, 'controller'), mode='w') as script:
            # shebang
            # if self.login_shell:
            #     script.write('#!/bin/bash -l\n')
            # else:
            #     script.write('#!/bin/bash\n')
            script.write('#!/usr/bin/env bash\n')

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
                script.write('export MPICH_RDMA_ENABLED_CUDA=1\n')
                script.write('export MPICH_G2G_PIPELINE=256\n')
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
                    if self.pgi_version == '16.9.0':
                        script.write('module unload pgi\n')
                        script.write('module load pgi/16.9.0\n')
                    elif self.pgi_version == '17.5.0':
                        script.write('module unload pgi\n')
                        script.write('module load pgi/17.5.0\n')
                    elif self.pgi_version == '17.10.0':
                        script.write('module unload pgi\n')
                        script.write('module use /apps/common/UES/pgi/17.10/modulefiles\n')
                        script.write('module load pgi/17.10\n')
                        script.write('export PGI_VERS_STR=17.10.0\n')
                    else:
                        raise ValueError("valid pgi versions are '16.9.0', '17.5.0' or '17.10.0'")

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
        ET.SubElement(config, 'machine').text = 'daint'
        ET.SubElement(config, 'account').text = self.account
        ET.SubElement(config, 'wall_time').text = self.wall_time
        ET.SubElement(config, 'partition').text = self.partition
        ET.SubElement(config, 'modules_opt').text = self.modules_opt
        ET.SubElement(config, 'pgi_version').text = self.pgi_version
        ET.SubElement(config, 'login_shell', type='py_eval').text = str(self.login_shell)


    def _submit_func(self):
        check_call(['sbatch', 'controller'])


    def _run_func(self):
        check_call(['module list'], shell=True)
        if self.cosmo_only:
            if self.gpu_mode:
                run_cmd = 'srun -u --ntasks-per-node=1 -n {:d} {:s}'.format(self._n_nodes, self.COSMO_exe)
                # check_call(['srun', '-u', '--ntasks-per-node=1', '-n', str(self._n_nodes), self.COSMO_exe])
            else:
                run_cmd = 'srun -u -n {:d} {:s}'.format(self._n_nodes * self._n_tasks_per_node, self.COSMO_exe)
                # check_call(['srun', '-u', '-n', str(self._n_nodes * self._n_tasks_per_node), self.COSMO_exe])
        else:
            self._build_proc_config()
            run_cmd = 'srun -u --multi-prog ./proc_config'
        print("running " + run_cmd)
        sys.stdout.flush()
        check_call(run_cmd, shell=True)

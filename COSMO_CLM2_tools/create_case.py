from __future__ import print_function
from .cc2_case import factory as cc2_case_factory, available_cases
from .tools import date_fmt, get_xml_node_args
from subprocess import check_call
from argparse import ArgumentParser, RawTextHelpFormatter, Action as arg_action
import f90nml
from datetime import datetime, timedelta
import os
import xml.etree.ElementTree as ET
import shutil

def create_case():
    """
    Create a Cosmo-CLM2 case from cmd line arguments and xml setup file

    See ``cc2_create_case --help``
    """

    # Build command line parser
    # =========================

    # Custom action factory to fill in cc2_cmd_args dictionnary
    cc2_cmd_args = {}
    case_actions = {}

    def cc2_act(*groups):

        for group in groups:
            if group not in cc2_cmd_args:
                cc2_cmd_args[group] = {}

        key = '.'.join(groups)

        if key not in case_actions:
            def call(self, parser, args, values, option_string=None):
                for group in self.cc2_groups:
                    cc2_cmd_args[group][self.dest] = values
            name = 'cc2_' + '_'.join(groups)
            case_actions[key] = type(name, (arg_action,),{'__call__': call, 'cc2_groups': groups})

        return case_actions[key]

    # function for boolean type
    def str_to_bool(val_str):
        return bool(eval(val_str))

    # Create parser
    dsc = "Set up and run a COSMO_CLM2 case\n"\
          "--------------------------------\n"\
          "Options can be set up either by xml file or the following command line arguments.\n"\
          "xml file options must be stored in a subelement of the root element tagged with 'main'.\n"\
          "and/or the specific machine (see https://github.com/COSMO-RESM/COSMO_CLM2_tools/blob/master/COSMO_CLM2_tools/default_setup.xml)\n"\
          "Command line arguments have precedence over xml file ones."
    parser = ArgumentParser(description=dsc, formatter_class=RawTextHelpFormatter)
    main_group = parser.add_argument_group('main', 'Options common to all machines')
    main_group.add_argument('--machine', metavar='MACH', action=cc2_act('main'),
                            help="machine on which the case is running (default: has to be given \n"\
                            "either by the command line or the xml setup file)")
    main_group.add_argument('-s', '--setup-file', metavar='FILE', help="xml file conatining setup options")
    main_group.add_argument('--name', action=cc2_act('main'), help="case name (default: COSMO_CLM2)")
    main_group.add_argument('--install_dir', action=cc2_act('main'),
                            help="directory where the case is installed (default: $SCRATCH on daint)")
    main_group.add_argument('--cosmo_only', action=cc2_act('main'), type=str_to_bool,
                            help="run only cosmo with build-in soil model TERRA\n"\
                            "(type: bool, using anything Python can parse as a boolean, default: False)\n"\
                            "Be carefull to provide a COSMO executable compiled accordingly")
    main_group.add_argument('--start_date', metavar='DATE_1', action=cc2_act('main'),
                            help="simulation start date formatted as YYYY-MM-DD-HH")
    main_group.add_argument('--end_date', metavar='DATE_2', action=cc2_act('main'),
                            help="simulation end date formatted as YYYY-MM-DD-HH")
    main_group.add_argument('--run_length', metavar='dt', action=cc2_act('main'),
                            help="sets simulation length if end_date not specified or run length\n"\
                            "between restarts otherwise\n"\
                            "dt is of the form 'N1yN2m' or 'N1y' or 'N2m' or 'N3d'\n"\
                            "N1, N2 and N3 being arbitrary integers (N2>12 possible) and\n"\
                            "'y', 'm' and 'd' stand for year, month and day")
    main_group.add_argument('--cos_in', action=cc2_act('main'),
                            help="COSMO input files directory (default: ./COSMO_input)")
    main_group.add_argument('--cos_nml', action=cc2_act('main'),
                            help="COSMO namelists directory (default: ./COSMO_nml)")
    main_group.add_argument('--cos_exe', action=cc2_act('main'),
                            help="path to COSMO executable (default: ./cosmo)")
    main_group.add_argument('--cesm_in', action=cc2_act('main'),
                            help="CESM input files directory (default: ./CESM_input)")
    main_group.add_argument('--cesm_nml', action=cc2_act('main'),
                            help="CESM namelists directory (default: ./CESM_nml)")
    main_group.add_argument('--cesm_exe', action=cc2_act('main'),
                            help="CESM executable (default: ./cesm.exe)")
    main_group.add_argument('--oas_in', action=cc2_act('main'),
                            help="OASIS input files directory (default: ./OASIS_input)")
    main_group.add_argument('--oas_nml', action=cc2_act('main'),
                            help="OASIS namelists directory (default: ./OASIS_nml)")
    main_group.add_argument('--ncosx', action=cc2_act('main'), type=int,
                            help="number of subdomains along the 'x-axis'\n"\
                            "for COSMO domain decomposition (type: int, default: from INPUT_ORG namelist)")
    main_group.add_argument('--ncosy', action=cc2_act('main'), type=int,
                            help="number of subdomains along the 'y-axis'\n"\
                            "for COSMO domain decomposition (type: int, default: from INPUT_ORG namelist)")
    main_group.add_argument('--ncosio', action=cc2_act('main'), type=int,
                            help="number of cores dedicated to i/o work'\n"\
                            "(type: int, default: from INPUT_ORG namelist)")
    main_group.add_argument('--ncesm', action=cc2_act('main'), type=int,
                            help="number of subdomains for CESM domain decomposition'\n"\
                            "(type: int, default: from drv_in namelist)")
    main_group.add_argument('--gpu_mode', action=cc2_act('main'), type=str_to_bool,
                            help="run COSMO on gpu (type: bool, using anything Python can parse as a boolean,\n"\
                            "default: False)")
    main_group.add_argument('--dummy_day', action=cc2_act('main'), type=str_to_bool,
                            help="perform a dummy day run after end of simulation to get last COSMO output.\n"\
                            "(type: bool, using anything Python can parse as a boolean, default: True)")

    slurm_group = parser.add_argument_group('slurm', 'Options specific to the slurm workload manager.\n'\
                                            '(common to all machines using the slurm scheduler)')
    slurm_group.add_argument('--wall_time', action=cc2_act('daint', 'mistral'),
                             help="reserved time on compute nodes\n"\
                             "(default: '24:00:00' on daint, '08:00:00' on mistral)")
    slurm_group.add_argument('--account', action=cc2_act('daint', 'mistral'),
                             help="account to use for batch script\n"\
                             "(default: infered from $PROJECT on daint, None on mistral)")
    slurm_group.add_argument('--partition', action=cc2_act('daint', 'mistral'),
                             help="select a queue (default: None)")

    daint_group = parser.add_argument_group('daint', 'Options specific to the Piz Daint machine')
    daint_group.add_argument('--modules_opt', action=cc2_act('daint'), choices=['switch', 'none', 'purge'],
                             help="Option for loading modules at run time (default: switch)")
    daint_group.add_argument('--pgi_version', action=cc2_act('daint'),
                             help="specify pgi compiler version at run time (default: None)")
    daint_group.add_argument('--shebang', action=cc2_act('daint'),
                             help="submit script shebang (default: #!/usr/bin/env bash)")

    cmd_line_group = parser.add_argument_group('cmd line', 'Options only avialble to the command line (no xml)')
    cmd_line_group.add_argument('--no_submit', action='store_false', dest='submit',
                                help="do not submit job after setup\n"\
                                "only command line argument, cannot be set in xml file")
    cmd_line_group.add_argument('--gen_oasis', action='store_true',
                                help="generate OASIS auxiliary files\n"\
                                "note that OASIS will crash after producing the files\n"\
                                "only command line argument, cannot be set in xml file\n")

    opts = parser.parse_args()

    # Parse machine and case argumennts from cmd line args and xml file
    # =================================================================
    machine, cc2_args = get_case_args(opts, cc2_cmd_args)
    print('- ML - DBG: cc2_args = ', cc2_args)

    # Create case instance
    # ====================
    cc2case = cc2_case_factory(machine, **cc2_args)

    # Change/delete namelists parameters following xml file
    # =====================================================
    modify_nml_from_xml(cc2case, opts)

    # Finalize
    # ========
    cc2case.write_open_nml()
    cc2case.to_xml(file_name='config.xml')

    # Submit case
    # ===========
    if opts.submit:
        cc2case.submit()

def get_case_args(cmd_opts, cc2_cmd_args):

    if cmd_opts.gen_oasis:
        cc2_cmd_args['main']['dummy_day'] = False

    if 'machine' in cc2_cmd_args['main']:
        machine = cc2_cmd_args['main']
    else:
        machine = None

    xml_file = cmd_opts.setup_file
    if xml_file:
        tree_root = ET.parse(xml_file).getroot()
        main_node = tree_root.find('main')
        if machine is None:
            machine_name_node = main_node.find('machine')
            if machine_name_node is not None:
                machine = machine_name_node.text
        machine_node = tree_root.find(machine)

    if machine is None:
        raise ValueError("'machine' option has to be given either by the command line or the xml setup file")

    main_args = get_xml_node_args(main_node, exclude=('machine'))
    main_args.update(cc2_cmd_args['main'])

    machine_args = get_xml_node_args(machine_node)
    machine_args.update(cc2_cmd_args[machine])

    cc2_args = {k:v for k,v in main_args.items() if v is not None}
    cc2_args.update({k:v for k,v in machine_args.items() if v is not None})
    cc2_args['install'] = True

    return machine, cc2_args

def modify_nml_from_xml(cc2case, cmd_opts):
    """Modify case namelists following instructions from xml setup file"""

    if cmd_opts.setup_file is None:
        return

    tree_root = ET.parse(cmd_opts.setup_file).getroot()

    # Change parameters
    nodes = tree_root.findall('change_par')
    if nodes:
        for node in nodes:
            name = node.get('file')
            block = node.get('block')
            n = node.get('n')
            param = node.get('param')
            val_str = node.text
            if name is None:
                raise ValueError("'file' xml attribute is required to change parameter")
            if block is None:
                raise ValueError("'block' xml attribute is required to change parameter")
            if param is None:
                raise ValueError("'param' xml attribute is required to change parameter")
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
                cc2case.nml[name][block][param] = value
            else:
                cc2case.nml[name][block][int(n)-1][param] = value

    # Delete parameters
    nodes = tree_root.findall('del_par')
    if nodes:
        for node in nodes:
            name = node.get('file')
            block = node.get('block')
            n = node.get('n')
            param = node.get('param')
            if name is None:
                raise ValueError("'file' xml attribute is required to delete parameter")
            if block is None:
                raise ValueError("'block' xml attribute is required to delete parameter")
            if param is None:
                raise ValueError("'param' xml attribute is required to delete parameter")
            if n is None:
                del cc2case.nml[name][block][param]
            else:
                del cc2case.nml[name][block][int(n)-1][param]

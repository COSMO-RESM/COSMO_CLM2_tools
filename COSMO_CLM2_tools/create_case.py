from __future__ import print_function
from .cc2_case import factory as cc2_case_factory, available_cases
from .date_formats import date_fmt_in, date_fmt_cosmo
from subprocess import check_call
from argparse import ArgumentParser, RawTextHelpFormatter
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

    # Parse setup options from command line and xml file
    # ==================================================

    # Options from command line
    # -------------------------
    def str_to_bool(val_str):
        return bool(eval(val_str))


    dsc = "Set up and run a COSMO_CLM2 case\n"\
          "--------------------------------\n"\
          "Options can be set up either by xml file or the following command line arguments.\n"\
          "xml file options must be stored in a subelement of the root element tagged with 'main'.\n"\
          "and/or the specific machine (see https://github.com/COSMO-RESM/COSMO_CLM2_tools/blob/master/COSMO_CLM2_tools/default_setup.xml)\n"\
          "Command line arguments have precedence over xml file ones."
    parser = ArgumentParser(description=dsc, formatter_class=RawTextHelpFormatter)
    main_group = parser.add_argument_group('main', 'Options common to all machines')
    main_group.add_argument('--machine', metavar='MACH', choices=['daint', 'mistral'],
                            help="machine on which the case is running (default: has to be given \n"\
                            "either by the command line or the xml setup file)")
    main_group.add_argument('-s', '--setup-file', metavar='FILE', help="xml file conatining setup options")
    main_group.add_argument('--name', help="case name (default: COSMO_CLM2)")
    main_group.add_argument('--path', help="directory where the case is set up (default: $SCRATCH/NAME on daint)")
    main_group.add_argument('--cosmo_only', type=str_to_bool, help="run only cosmo with build-in soil model TERRA\n"\
                            "(type: bool, using anything Python can parse as a boolean, default: False)\n"\
                            "Be carefull to provide a COSMO executable compiled accordingly")
    main_group.add_argument('--start_date', metavar='DATE_1',
                            help="simulation start date formatted as YYYY-MM-DD-HH")
    main_group.add_argument('--end_date', metavar='DATE_2',
                            help="simulation end date formatted as YYYY-MM-DD-HH")
    main_group.add_argument('--run_length', metavar='dt',
                            help="sets simulation length if end_date not specified or run length\n"\
                            "between restarts otherwise\n"\
                            "dt is of the form 'N1yN2m' or 'N1y' or 'N2m' or 'N3d'\n"\
                            "N1, N2 and N3 being arbitrary integers (N2>12 possible) and\n"\
                            "'y', 'm' and 'd' stand for year, month and day")
    main_group.add_argument('--cos_in', help="COSMO input files directory (default: ./COSMO_input)")
    main_group.add_argument('--cos_nml', help="COSMO namelists directory (default: ./COSMO_nml)")
    main_group.add_argument('--cos_exe', help="path to COSMO executable (default: ./cosmo)")
    main_group.add_argument('--cesm_in', help="CESM input files directory (default: ./CESM_input)")
    main_group.add_argument('--cesm_nml', help="CESM namelists directory (default: ./CESM_nml)")
    main_group.add_argument('--cesm_exe', help="CESM executable (default: ./cesm.exe)")
    main_group.add_argument('--oas_in', help="OASIS input files directory (default: ./OASIS_input)")
    main_group.add_argument('--oas_nml', help="OASIS namelists directory (default: ./OASIS_nml)")
    main_group.add_argument('--ncosx', type=int, help="number of subdomains along the 'x-axis'\n"\
                            "for COSMO domain decomposition (type: int, default: from INPUT_ORG namelist)")
    main_group.add_argument('--ncosy', type=int, help="number of subdomains along the 'y-axis'\n"\
                            "for COSMO domain decomposition (type: int, default: from INPUT_ORG namelist)")
    main_group.add_argument('--ncosio', type=int, help="number of cores dedicated to i/o work'\n"\
                            "(type: int, default: from INPUT_ORG namelist)")
    main_group.add_argument('--ncesm', type=int, help="number of subdomains for CESM domain decomposition'\n"\
                            "(type: int, default: from drv_in namelist)")
    main_group.add_argument('--gpu_mode', type=str_to_bool,
                            help="run COSMO on gpu (type: bool, using anything Python can parse as a boolean,\n"\
                            "default: False)")
    main_group.add_argument('--dummy_day', type=str_to_bool,
                            help="perform a dummy day run after end of simulation to get last COSMO output.\n"\
                            "(type: bool, using anything Python can parse as a boolean, default: True)")

    slurm_group = parser.add_argument_group('slurm', 'Options specific to the slurm workload manager')
    slurm_group.add_argument('--wall_time', help="reserved time on compute nodes\n"\
                             "(default: '24:00:00' on daint, '08:00:00' on mistral)")
    slurm_group.add_argument('--account', help="account to use for batch script\n"\
                             "(default: infered from $PROJECT on daint, None on mistral)")
    slurm_group.add_argument('--partition', help="select a queue (default: None)")

    daint_group = parser.add_argument_group('daint', 'Options specific to the Piz Daint machine')
    daint_group.add_argument('--modules_opt', choices=['switch', 'none', 'purge'],
                             help="Option for loading modules at run time (default: switch)")
    daint_group.add_argument('--pgi_version', help="specify pgi compiler version at run time (default: None)")
    daint_group.add_argument('--shebang', help="submit script shebang (default: #!/usr/bin/env bash)")

    cmd_line_group = parser.add_argument_group('cmd line', 'Options only avialble to the command line (no xml)')
    cmd_line_group.add_argument('--no_submit', action='store_false', dest='submit',
                                help="do not submit job after setup\n"\
                                "only command line argument, cannot be set in xml file")
    cmd_line_group.add_argument('--gen_oasis', action='store_true',
                                help="generate OASIS auxiliary files\n"\
                                "note that OASIS will crash after producing the files\n"\
                                "only command line argument, cannot be set in xml file\n")

    opts = parser.parse_args()
    if opts.gen_oasis:
        opts.dummy_day = False

    # Set options to xml value if needed or default if nothing provided then perform some checks
    # ------------------------------------------------------------------------------------------
    # options defaults
    defaults = {'main': {'machine': None, 'name': 'COSMO_CLM2', 'path': None,
                         'cosmo_only': False, 'gen_oasis': False,
                         'start_date': None, 'end_date': None, 'run_length': None,
                         'cos_in': './COSMO_input', 'cos_nml': './COSMO_nml', 'cos_exe': './cosmo',
                         'cesm_in': './CESM_input', 'cesm_nml': './CESM_nml', 'cesm_exe': './cesm.exe',
                         'oas_in': './OASIS_input', 'oas_nml': './OASIS_nml',
                         'ncosx': None, 'ncosy': None, 'ncosio': None, 'ncesm': None,
                         'dummy_day': True, 'gpu_mode': False},
                'daint': {'wall_time': '24:00:00', 'account': None, 'partition': None,
                          'modules_opt': 'switch', 'pgi_version': None, 'shebang': '#!/usr/bin/env bash'},
                'mistral': {'wall_time': '08:00:00', 'account': None, 'partition': None}}

    # Apply default main options
    if opts.setup_file is not None:
        tree = ET.parse(opts.setup_file)
        xml_node = tree.getroot().find('main')
    else:
        xml_node = None
    apply_defaults(opts, xml_node, defaults['main'])

    # Check machine
    if opts.machine is None:
        raise ValueError("'machine' option has to be given either by the command line or the xml setup file")
    elif opts.machine not in available_cases:
        raise ValueError("invalid 'machine' option. Has to be either of " + str(list(available_cases.keys())))

    # Apply default machine specific options
    if opts.machine not in defaults:
        raise NotImplementedError("default options not implemented for machine {:s}".format(opts.machine))
    else:
        if opts.setup_file is not None:
            tree = ET.parse(opts.setup_file)
            xml_node = tree.getroot().find(opts.machine)
        else:
            xml_node = None
        apply_defaults(opts, xml_node, defaults[opts.machine])

    # Check case path
    # - ML - that will go to the daint_case class as soon as cleaner defaults are implemented
    if opts.path is None:
        if opts.machine == 'daint':
            opts.path = os.path.join(os.environ['SCRATCH'], opts.name)
        else:
            raise NotImplementedError("default path not implemented for machine {:s}".format(opts.machine))

    # Log
    # ===
    log = 'Setting up case {:s} in {:s}'.format(opts.name, opts.path)
    print(log + '\n' + '-' * len(log))

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
    if not opts.cosmo_only:
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
    # base arguments
    case_args = {'name': opts.name, 'path': opts.path,
                 'cosmo_only': opts.cosmo_only, 'gen_oasis': opts.gen_oasis,
                 'start_date': opts.start_date, 'end_date': opts.end_date,
                 'run_length': opts.run_length, 'COSMO_exe': os.path.basename(opts.cos_exe),
                 'CESM_exe': os.path.basename(opts.cesm_exe), 'ncosx': opts.ncosx,
                 'ncosy': opts.ncosy, 'ncosio': opts.ncosio, 'ncesm': opts.ncesm,
                 'gpu_mode': opts.gpu_mode, 'dummy_day': opts.dummy_day}
    # machine specific arguments
    if opts.machine == 'daint':
        machine_args = {'wall_time': opts.wall_time, 'account': opts.account,
                        'partition': opts.partition,'modules_opt': opts.modules_opt,
                        'pgi_version': opts.pgi_version, 'shebang': opts.shebang}
    elif opts.machine == 'mistral':
        machine_args = {'wall_time': opts.wall_time, 'account': opts.account,
                        'partition': opts.partition}
    else:
        raise NotImplementedError("machine_args dict not implemented for machine {:s}".format(opts.machine))
    case_args.update(machine_args)
    # create case instance
    cc2case = cc2_case_factory(opts.machine, **case_args)

    # Change/delete namelists parameters following xml file
    # =====================================================
    if opts.setup_file is not None:
        # Modify namelist parameters
        nodes = tree.getroot().findall('change_par')
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

        # Delete namelist parameters
        nodes = tree.getroot().findall('del_par')
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

    # Finalize
    # ========
    cc2case.write_open_nml()
    cc2case.to_xml(file_name='config.xml')

    # Submit case
    # ===========
    if opts.submit:
        cc2case.submit()

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
                    if opt_val_str is None:
                        apply_def = True
                    else:
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
                                                 + " is not a valid python type nor 'py_eval'")
        if apply_def:
            setattr(opts, opt, default)

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

    if os.path.exists('transfer_list'):
        os.remove('transfer_list')

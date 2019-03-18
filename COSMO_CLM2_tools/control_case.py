from .cc2_case import factory as cc2_case_factory
from .tools import get_xml_node_args
from argparse import ArgumentParser, RawTextHelpFormatter
import xml.etree.ElementTree as ET
from time import sleep


def control_case():
    # Parse arguments
    dsc = "Control a COSMO_CLM2 case"
    parser = ArgumentParser(description=dsc, formatter_class=RawTextHelpFormatter)
    parser.add_argument('xml_path', help="path to xml file containing case description")
    parser.add_argument('--action', choices=['run', 'transfer'], default='run',
                        help="path to xml file containing case description")
    cfg = parser.parse_args()

    # build cc2case object from xml file
    config = ET.parse(cfg.xml_path).getroot()
    machine = config.find('machine').text
    case_args = get_xml_node_args(config.find('main'))
    case_args.update(get_xml_node_args(config.find(machine)))
    cc2case = cc2_case_factory(machine, **case_args)

    if cfg.action == 'run':
        cc2case.run_status = 'running'

        # Submit next transfer
        if (cc2case._run_end_date < cc2case.end_date and cc2case.transfer_by_chunck):
            cc2case.transfer_status = 'submitted'
            cc2case.submit_next_transfer()

        # Run
        cc2case.run()
        cc2case.set_next_run()

        # Archive
        if cc2case.archive_dir is not None:
            cc2case.submit_archive()

        cc2case.run_status = 'complete'

        # Submit next run
        if (cc2case._run_end_date < cc2case.end_date and cc2case.transfer_status == 'complete'):
            cc2case.run_status = 'submitted'
            cc2case.submit_next_run()

    elif cfg.action == 'transfer':
        # Transfer
        cc2case.transfer_status = 'transferring'
        cc2case.transfer_input()
        cc2case.transfer_status = 'complete'

        # Submit next run
        if cc2case.run_status == 'complete':
            # avoid (limit a lot) race conditions
            sleep(10)
            if cc2case.run_status == 'complete':
                cc2case.run_status = 'submitted'
                cc2case.submit_next_run()

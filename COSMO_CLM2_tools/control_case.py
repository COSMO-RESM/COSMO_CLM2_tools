from .cc2_case import factory as cc2_case_factory
from .tools import get_xml_node_args
from argparse import ArgumentParser, RawTextHelpFormatter
import xml.etree.ElementTree as ET


def control_case():
    # Parse arguments
    dsc = "Control a COSMO_CLM2 case"
    parser = ArgumentParser(description=dsc, formatter_class=RawTextHelpFormatter)
    parser.add_argument('xml_path', help="path to xml file containing case description")
    cfg = parser.parse_args()

    # build cc2case object from xml file
    config = ET.parse(cfg.xml_path).getroot()
    machine = config.find('machine').text
    case_args = get_xml_node_args(config.find('main'))
    case_args.update(get_xml_node_args(config.find(machine)))
    cc2case = cc2_case_factory(machine, **case_args)

    # Run
    cc2case.run()

    # Submit next run
    if cc2case.set_next_run():
        cc2case.submit()

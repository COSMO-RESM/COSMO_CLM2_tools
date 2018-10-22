from .base_case import base_case
from argparse import ArgumentParser, RawTextHelpFormatter
import os


def control_case():
    # Parse arguments
    dsc = "Control a COSMO_CLM2 case"
    parser = ArgumentParser(description=dsc, formatter_class=RawTextHelpFormatter)
    parser.add_argument('xml_path', help="path to xml file containing case description")
    cfg = parser.parse_args()

    # Read case configuration from xml file
    path, xml_file = os.path.split(cfg.xml_path)
    os.chdir(path)
    cc2case = base_case.from_xml(xml_file)

    # Run
    cc2case.run()

    # Submit next run
    if cc2case.set_next_run():
        cc2case.submit()

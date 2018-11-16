from .cc2_case import factory as cc2_case_factory
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
                                     + " has to be a valid python type or 'py_eval'")


    cc2case = cc2_case_factory(machine, **args)

    # Run
    cc2case.run()

    # Submit next run
    if cc2case.set_next_run():
        cc2case.submit()

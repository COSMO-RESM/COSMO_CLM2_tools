date_fmt = {'in': '%Y-%m-%d-%H', 'cosmo': '%Y%m%d%H','cesm': '%Y%m%d'}


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


def get_xml_node_args(node, exclude=()):
    """Read case arguments from xml node"""

    if node is None:
        return {}

    xml_args = {}

    for opt in node.iter():
        if opt is not node and opt.tag not in exclude:
            if opt.get('type') is None:
                xml_args[opt.tag] = opt.text
            elif opt.get('type') == 'py_eval':
                xml_args[opt.tag] = eval(opt.text)
            else:
                opt_type = eval(opt.get('type'))
                if isinstance(opt_type, type):
                    xml_args[opt.tag] = opt_type(opt.text)
                else:
                    raise ValueError("xml atribute 'type' " + opt.get('type')
                                     + " for node " + opt.tag
                                     + " has to be a valid python type or 'py_eval'")

    return xml_args

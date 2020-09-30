import re


# This function turns user input into x,y coordinates,
# even if the user is not careful about syntax.
def parse_user_xy(user_xy):
    # The user may have included extra spaces.
    xy = str(user_xy).strip()
    xy = re.sub(r'\s+', ' ', xy)

    # The user may have omitted commas in lists.
    xy_string = ''
    for st in range(len(xy)):
        if st == 0 or st == len(xy) - 1:
            xy_string += xy[st]
        elif xy[st] == ' ' and xy_string[-1].isdigit():
            if xy[st + 1].isdigit() or xy[st + 1] == '.' or xy[st + 1] == '-':
                xy_string += ','
        else:
            xy_string += xy[st]

    # The user may have omitted parentheses.
    parenths = False
    for st in xy:
        if st == '(' or st == '[':
            parenths == True
            break
    if parenths == False:
        xy = ('(' + xy_string + ')')
    else:
        xy = xy_string

    # Evaluate the user string, forcing a list even if it evaluates as tuple.
    return list(eval(xy))

def parse_cct(cct):
    """
    Parses a user string representing a correlated color temperature (CCT)
    (possibly including a trailing 'K') into a float

    Parameters
    ----------
    cct_ : str
        Correlated color temperature as a string of some sort

    Returns
    -------
    float

    Notes
    -----
    Differs from Steve's parser in that it won't ignore random crap in the
    middle of a number, won't allow multiple decimal points, won't allow
    CCTs below 1 (e.g. 0.341), and doesn't allow negative CCTs.

    That leaves: 5400, 5400., 5400.1, 5400K, 5400.K, 5400.1K, 5600 K, and the like

    TODO: see if original parser would have returned a complex number for, e.g., 5400-2i
    """
    p = re.compile(r'^\s*(\d+(\.\d*)?)\s*K?\s*$')
    m = p.match(cct)
    if m:
        return float(m.group(1))
    else:
        raise SyntaxError(f"could not parse {cct} as a correlated color temperature")

def parse_tint(tint):
    """
    Parses a user string representing a tint. Note that tint is expressed here in
    delta uv units and not in the units found on (say) an ALEXA or in the units
    proposed for the ACES container spec (SMPTE ST 2065-4:2021)

    Parameters
    ----------
    tint : str
        tint expressed in delta uv

    Returns
    -------
    float

    Notes
    -----
    Differs from Steve's parser in that it won't ignore random crap in the
    middle of a number, won't allow multiple decimal points, won't allow
    CCTs below 1 (e.g. 0.341), and doesn't allow negative CCTs.

    That leaves: 5400, 5400., 5400.1, 5400K, 5400.K, 5400.1K, 5600 K, and the like

    TODO: see if original parser would have returned a complex number for, e.g., 5400-2i
    """
    p = re.compile(r'^\s*([-+]?\d+(\.\d*)?)\s*$')
    m = p.match(tint)
    if m:
        return float(m.group(1))
    else:
        raise SyntaxError(f"could not parse {tint} as a tint (expressed in delta uv units")


# Convert the user's white balance input to usable CCT and Tint.
def parse_user_cct_and_tint(user_cct, user_tint):
    if type(user_cct) == str:
        cct = ''
        for c in user_cct:
            if c.isdigit() or c == '.' or c == '-':
                cct += c
        cct = eval(cct)
    else:
        cct = float(user_cct)

    if type(user_tint) == str:
        tint = ''
        for c in user_tint:
            if c.isdigit() or c == '.' or c == '-':
                tint += c
        tint = eval(tint)
    else:
        tint = float(user_tint)

    return [cct, tint]

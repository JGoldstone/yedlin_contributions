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

    """
    p = re.compile(r'^\s*([-+]?\d+(\.\d*)?)\s*$')
    m = p.match(tint)
    if m:
        return float(m.group(1))
    else:
        raise SyntaxError(f"could not parse {tint} as a tint (expressed in delta uv units")


# Convert the user's white balance input to usable CCT and Tint.
def parse_user_cct_and_tint(user_cct, user_tint):
    return [parse_cct(user_cct), parse_tint(user_tint)]


def user_message():
    print('\n\n-----------------------------------'
          '\nACES Data Set Calculator v0.01\n'
          'This calculator returns the expected ACES rgb triplets for a Macbeth chart\n'
          'photographed with various illuminants and with various user settings for\n'
          'white balance.\n'
          '\n-----------------------------------\n'
          '\nFor guidance on how to shoot a data set, enter \'help\' or \'h\'\n'
          'Or enter anything else to continue.\n')
    user_help = input()

    uh = user_help.lower()
    if uh == 'help' or uh == 'h':
        print('''
RECOMMENDATIONS FOR SHOOTING AN ACES IDT DATA SET:

MATERIALS:

-The camera you are testing.

-A high quality, contemporary lens. Medium to long focal length. Recently
serviced, including verification of t/stop markings\' accuracy. Use a prime lens,
not a zoom.

-A Macbeth Color Checker chart. In good condition, not faded.

-A Kodak R-27 18% Gray Card. In good condition, not faded.

-An exposure meter that can take reflective spot readings, recently calibrated
(using calibration constant K = 12.5)

-A spectrometer or spectroradiometer, calibrated. (Not a colorimeter.)

-2 luminaires (lighting instruments). Both luminaires should be of the same
model and in the same condition. They should be as close to the whitepoints you 
want to test as possible. Or they can be filtered/gelled to close to that 
whitepoint. They should be smooth spectrum illuminants such as incandescents,
not LEDs.

-Support equipment as needed (tripod, light stands, etc.)

-A room that can be light controlled, so when the work lights are off there is
no light contamination -- make sure your luminaires are the only lights in the
room when taking readings or rolling the camera.


METHOD:

1.   Arrange the camera so that it is shooting the chart square on and so that
the chart is in the center of the frame. Adjust the camera\'s distance from the
chart so that the chart takes up about one quarter to one third of the frame's
horizontal. Do not get closer, in order to avoid vignetting of the lens. If the
lens's image circle does not fully cover the camera's imager, then the chart
should be even smaller in the frame: a third to a quarter of the image circle's
diameter.

2.   Make sure that the Macbeth chart is flat and not curved or warped. Place
the R-27 gray card on top of it (you'll want to be able to meter the gray card 
in the hero postion, but be able to remove it and reveal the Macbeth chart
without disturbing the assembly). Make sure that the gray card is also flat 
and not curved or warped.

3.   Arrange the luminaires on either side of the camera, at least 45 degrees off
of the lens axis, so that they are symmetrical and evenly illuminating the chart 
Make sure there are no sheens or caustic reflections from the chart to the
line-of-sight to the camera. The camera should only receive diffused light from
the chart -- no glares.

4.   Set the camera and the exposure meter to matching ISO, shutter angle, and
frame-rate settings.

5.   Stand in the line-of-sight of the camera and use the exposure meter to take 
a spot reading of the R-27 gray card. It's important to use the R-27 card and 
not the Macbeth chart in this step, becuase the Macbeth\'s chart "middle" gray
is 19.8% reflectance, not 18%.

6.   Based on the reading, adjust the luminaires, the camera\'s shutter angle,
and/or the camera's frame rate so that the f/stop reading from the exposure meter
is precisely on an f/stop that\'s marked on the lens (for example, f/4.0 but not
f/4.3). Also make sure that the reading is an aperture no wider than f/4.0 and
at least 1 stop closed down from wide open on your lens. If you change the shutter
or frame rate on the camera, match the changes in the exposure meter\'s settings.*

7.   Keep taking new readings and repeating step 6 until all criteria are met,
taking care to meter all over the gray card, making sure that it's evenly
illuminated.

8.   Set the lens\'s t/stop ring to the f/stop reading from the exposure meter.
When setting the stop, start with the iris wide open and close down to the marked
t/stop -- don't start closed and open to your stop (becuase lenses are marked
to account for backlash in the mechanism). Doublecheck that the ISO, shutter angle, 
and frame rate settings in the camera are the same as those in the meter and that
the f/stop spot-reading of the gray card matches the iris\'s t/stop setting.
Also doublecheck again that the gray card is evenly illuminated all over.

9.   Remove the R-27 card to reveal the Macbeth chart and roll the camera to
create your data set for analysis.

10.   Remove the lens from the camera, and set the aperture ring to one stop closed
from wide open. Set your spectrometer or spectroradiometer to take an emissive
reading. Hold the lens between the spectrometer and one of the luminaires and take
a reading, noting the CIE 1931 x,y chromaticity coordinates. Repeat for the other
luminaire. Compute the average x,y coordinate of the two readings. This yields the
data set's white point as seen by the camera through the lens.

* This recommended method includes shooting only one neutral chart for
creating or verifying chromatic adaptation of an IDT for one specific ISO. However, 
in reality, you will probably be gathering data for multiple ISO settings and data
for power-transfer functions at the same time. So, you will likely also want to use
the shutter and framerate to shoot an exposure sweep of over- and under- exposure
around this neutral chart. Correspondingly, in step 6, it is adviasable to adjust
the luminaires and the f/stop so that the desired exposure is achieved in the middle
of the camera's frame-rate and shutter-angle ranges, for maximum latitude to adjust
both directions up and down. For most cameras at the time of writing, this will be
soemwhere around 16 frames per second and 45-degree shutter.
        ''')

def user_interface():
    print('-----------------------------------\n\n\n'
          'Enter the CIE 1931 x,y coordinates of the illuminant that fell on the Macbeth chart.\n'
          'Mesure the illumant through the same lens you used to shoot the chart, so that\n'
          'any color shift from the lens can be cancelled out in this calculation.\n\n'
          'x,y coordinant of illuminant, through the lens?')
    illuminant = input()

    print('\nEnter the camera white balance CCT in degrees Kelvin that you are testing\n'
          'or for which you\'re trying to build an IDT.\n'
          'User CCT setting?')
    cct = parse_cct(input())

    print('\nEnter the camera white balnce tint in ', end='')
    try:
        print(delta, end='')
    except:
        print('Delta-', end='')
    print('uv that you are testing or for which you\'re trying to build an IDT.\n'
          'User tint setting?')
    tint = parse_tint(input())

    print('\nDo you want the targets for the whole Macbeth chart or just for\n'
          'middle gray, red, green, blue, cyan, magenta, and yellow?\n'
          '\t1. 18% gray and Whole Macbeth chart\n'
          '\t2. 18% gray and Macbeth Gray + RGBCMY')
    chips = input()
    chips = chips.strip()

    print('\nEnter "n" if you would like results as a Nuke node.\n'
          'otherwise, results will be given as human-readable list.')

    nuke = input()
    nuke = nuke.strip()
    if nuke.lower() == 'n':
        nuke = True
    else:
        nuke = False

    if chips[0] == '2':
        colorChips = [macbeth[0]] + [macbeth[22]] + macbeth[13:19]
    else:
        colorChips = macbeth

    print('\nIn what format would you like your ACES values returned?\n'
          '    1. ACES AP0 linear\n'
          '    2. ACES AP1 linear\n'
          '    3. ACEScc log\n'
          '    4. ACEScct log\n')
    ACES_format = input()
    ACES_format = ACES_format.strip()

    if ACES_format[0] == '2':
        colorspace = 'AP1'
    elif ACES_format[0] == '3':
        colorspace = 'CC'
    elif ACES_format[0] == '4':
        colorspace = 'CCT'
    else:
        colorspace = 'AP0'
    return (illuminant, cct, tint, colorChips, colorspace)

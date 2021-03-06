#!/usr/bin/python
from __future__ import print_function

#ACES Data Set Calculator
#by Steve Yedlin, ASC
#First Version: v0.01 March 22, 2018
#This Version: v0.02 April 25, 2018
#Fixed mistake in "18% gray" line on Sept 1, 2020

import re
import numpy as np
from collections import namedtuple
from copy import deepcopy
from math import log

try:
  input = raw_input
except:
  pass

# First, define some component color conversions.
# Not all of these are in use in the current iteration of this script,
# but are here for future use. These have been cross-referenced to find agreement
# with Giorgianni/Madden, brucelindbloom.com, Nuke, easyrgb.com, Colour Science 
# For Python, and Wikipedia

def Yxy_to_XYZ(triplet):

    Y,x,y = triplet
    XYZ = []

    if y == 0: XYZ.append(0.0)
    else : XYZ.append(x * ( Y / y))

    XYZ.append(Y)

    if y == 0: XYZ.append(0.0)
    else: XYZ.append(( 1 - x - y ) * ( Y / y ))

    return XYZ


def XYZ_to_Yxy(triplet):

    X,Y,Z = triplet
    Yxy = []

    Yxy.append(Y)

    if(X+Y+Z) == 0: Yxy.append(0.0)
    else : Yxy.append(X / (X + Y + Z))

    if (X+Y+Z) == 0: Yxy.append(0.0)
    else: Yxy.append (Y / ( X + Y + Z ))

    return Yxy
    
#Matrix coefficients per Edward Giorgianni and Thomas Madden in "Digital Color Management"
#match Lindbloom's to at least 3 decimal places.
def XYZ_to_linRGB(triplet):
    X,Y,Z = triplet
    R = X*3.242 + Y*-1.5374 + Z*-.4986
    G = X*-.9692 + Y*1.876 + Z*.0416
    B = X*.0556 + Y*-.2040 + Z*1.0570
    return [R,G,B]
    
#Inverse matrix derived from above published matrix
#matches Lindbloom to at least 3 decimal places.
def linRGB_to_XYZ(triplet):
    R,G,B = triplet
    X = R*0.4122109 + G*0.35742544 + B*0.18037791
    Y = R*0.21253222 + G*0.71506279 + B*0.07211159
    Z = R*0.01933552 + G*0.11920526 + B*0.95050308
    return [X,Y,Z]


def linRGB_to_sRGB(triplet):

    sRGB = []
    for val in triplet:
        if val > 0.0031308: sRGB.append(1.055 * (val ** ( 1 / 2.4 )) - 0.055)
        else: sRGB.append(12.92 * val)
    return sRGB

def sRGB_to_linRGB(triplet):

    linRGB = []
    for val in triplet:
        if val > 0.04045: linRGB.append((( val + 0.055 ) / 1.055 ) ** 2.4 )
        else: linRGB.append(val/12.92)
    return linRGB


def xy_to_uv(coords):
    x,y = coords

    u = 4*x / (12*y - 2*x + 3)
    v = 6*y / (12*y - 2*x + 3)

    return [u,v]

def uv_to_xy(coords):
    u,v = coords

    x = 3*u / (2*u - 8*v + 4)
    y = 2*v / (2*u - 8*v + 4)

    return [x,y]


# Here are some conversions between x,y and Kelvin (CCT) 
# that only work if the color is on the Planckian locus
def k_to_xy (k):

    if k <= 4000:
        x = -.2661239*(10**9)/(k**3)-.2343580*(10**6)/(k**2)+.8776956*(10**3)/(k)+.179910
    else: 
        x = -3.0258469*(10**9)/(k**3)+2.1070379*(10**6)/(k**2)+.2226347*(10**3)/(k)+.240390

    if k <= 2222:
        y = -1.1063814*(x**3)-1.34811020*(x**2)+2.18555832*x-.20219683
    elif k <= 4000:
        y = -.9549476*(x**3)-1.37418593*(x**2)+2.09137015*x-.16748867
    else:
        y = 3.0817580*(x**3)-5.87338670*(x**2)+3.75112997*x-.37001483


    return [x,y]

def xy_to_k (xy):

    xVals = []
    kVals = []

    for k in range(1667 , 2060, 20):
        xVals.append(k_to_xy(k)[0])
        kVals.append(float(k))
    for k in range(2100 , 250010 , 100):
        xVals.append(k_to_xy(k)[0])
        kVals.append(float(k))

    xVals.sort()
    kVals.sort(reverse = True)

    k = np.interp(xy[0] , xVals , kVals) 

    if xy[0] > 0.56463 or xy[0] < 0.2524:
        k = 0.0

    return k

# CAT02 is the hero chomatic adaptation method, not Bradford, but here's Bradford anyway.
# Bradford Matrix coefficients as given by Bruce Lindblom at http://www.brucelindbloom.com
def XYZ_to_LMS(triplet):
    X,Y,Z = triplet
    L = X*0.8951000 + Y*0.2664000 + Z*-0.1614000
    M = X*-0.7502000 + Y*1.7135000 + Z*0.0367000
    S = X*0.0389000 +Y*-0.0685000 + Z*1.0296
    return [L,M,S]


#These inverse Bradford coefficients are derived, not from published documents.
def LMS_to_XYZ(triplet):
    X,Y,Z = triplet
    L = X*0.986993 + Y*-0.147054 + Z*0.159963
    M = X*0.432305 + Y*0.51836 + Z*0.0492912
    S = X*-0.00852866 + Y*0.0400428 + Z*0.968487
    return [L,M,S]
    
# Convert XYZ to/from CAT02 using the CIECAM02 method.
# The coefficients are from Fairchild 'A Revision  of CIECAM97 for Practical Application'
# at: http://rit-mcsl.org/fairchild/PDFs/AppearanceLec.pdf
# Wikipedia quotes this document and other Fairchild work in CAT02 entries.
def XYZ_to_CAT02(triplet):
    X,Y,Z = triplet
    L = X*0.7328 + Y*0.4296 +Z*-0.1624
    M = X*-0.7036 + Y*1.6975 + Z*0.0061
    S = X*0.0030 + Y*0.0136 +Z*0.9834
    return [L,M,S]
    
#This is a derived inverse for more decimal places than Fairchild's published inverse.
def CAT02_to_XYZ(triplet):
    L,M,S = triplet
    X = L*1.09612 + M*-0.278869 + S*0.182745
    Y = L*0.454369 + M*0.473533 + S*0.0720978
    Z = L*-0.00962761 + M*-0.00569803 + S*1.01533
    return [X,Y,Z]

# These conversions are for CIE L*a*b*, not the pre-CIE Lab,
# but we're not using asterisks in the function name.
def XYZ_to_Lab (triplet, refWhite=[.95045,1.,1.08892]):

    refX, refY, refZ = refWhite
    
    x = triplet[0]/refX
    y = triplet[1]/refY
    z = triplet[2]/refZ
    
    if x > 0.008856: x = x ** (1.0/3.0)
    else: x = (7.787 * x) + (16.0 / 116.0)
    
    if (y > 0.008856): y = y ** (1.0/3.0)
    else: y = (7.787 * y) + (16.0 / 116.0)
    
    if (z > 0.008856): z = z ** (1.0/3.0)
    else: z = (7.787 * z) + (16.0 / 116.0)
    
    L = (116.0 * y) - 16.0
    a = 500.0 * (x - y)
    b = 200.0 * (y - z)
    
    #Return a normalized version.    
    return [L*.01,a*.01,b*.01]
    
def Lab_to_XYZ (triplet, refWhite=[.95045,1.,1.08892]):

    refX, refY, refZ = refWhite
    
    #We're using normalized CIE L*a*b*, so we need to de-normalize.
    L = triplet[0]*100.
    a = triplet[1]*100.
    b = triplet[2]*100.
    
    y = (L + 16.0) / 116.0
    x = a / 500.0 + y
    z = y - b / 200.0
    
    if (y**3  > 0.008856): y = y**3.0
    else: y = (y - 16.0 / 116.0) / 7.787
    
    if (x**3  > 0.008856): x = x**3.0
    else: x = (x - 16.0 / 116.0) / 7.787
    
    if (z**3  > 0.008856): z = z**3.0
    else: z = (z - 16.0 / 116.0) / 7.787
        
    x *= refX
    y *= refY
    z *= refZ
    
    return [x,y,z]

def xy_to_CAT02(coords):
    Yxy = [1.] + coords
    XYZ = Yxy_to_XYZ(Yxy)
    return XYZ_to_CAT02(XYZ)

def CAT02_to_Yxy(triplet):
    XYZ = CAT02_to_XYZ(triplet)
    return XYZ_to_Yxy(triplet)
    
#The next group of conversions are based on the ACES documentation. 
def XYZ_to_linAP0(triplet):
    X,Y,Z= triplet
    r = (1.0498110175 * X) + (-0.0000974845*Z)
    g = (-0.4959030231*X) + (1.3733130458*Y) + (0.0982400361*Z)
    b = 0.9912520182 * Z
    return [r,g,b]
    
def linAP0_to_XYZ(triplet):
    r,g,b = triplet
    X = (r*0.9525523959) +  (b*0.0000936786)
    Y = (r*0.3439664498) + (g*0.7281660966) + (b*-0.0721325464)
    Z = b*1.0088251844
    return [X,Y,Z]
    
def linAP0_to_linAP1(triplet):
    r,g,b = triplet
    R = (r*0.6954522414) + (g*0.1406786965) + (b*0.16386906223)
    G = (r*0.0447945634) + (g* 0.8596711185) +(b*0.095534318275)
    B = (r*-0.0055258826) + (g*0.0040252103) + (b*1.0015006723)
    return [R,G,B]
    
def linAP1_to_linAP0(triplet):
    r,g,b = triplet
    R = (r*1.4514393161) + (g*-0.2365107469) + (b*-0.214928569337)
    G = (r*-0.0765537734) + (g*1.1762296998) + (b*-0.099675926475)
    B = (r*-0.0083161484) + (g*-0.0060324498) + (b*0.9977163014)
    return [R,G,B]
    
def linAP1_to_ACEScc(triplet):
    ACEScc = []
    for channel in triplet:
        if channel <= 0:
            ACEScc.append((log(2.**-.16,2.) + 9.72) / 17.52)
        elif channel < 2.**(-15.):
            ACEScc.append((log( (2.**-.16) + channel*.5 , 2.) + 9.72) / 17.52)
        else:
            ACEScc.append((log(channel,2.) + 9.72 ) / 17.52)
    return ACEScc
    
def ACEScc_to_linAP1(triplet):
    AP1 = []
    for channel in triplet:
        if channel <= (9.72-15)/17.52: 
            c = (2 ** ( channel*17.52 - 9.52 )) - (2.0**-16.)
            AP1.append( c*2.)
        elif channel < (log(65504.,2.) + 9.72) / 17.52:
            AP1.append (2. ** (channel*17.52 - 9.72))
        else: AP1.append(65504.)
    return AP1 
    
def linAP1_to_ACEScct(triplet):
    ACEScct = []
    for channel in triplet:
        if channel <=  0.0078125:
            ACEScct.append(10.5402377416545 * channel + 0.0729055341958355)
        else:
            ACEScct.append((log(channel,2.) + 9.72) / 17.52)
    return(ACEScct)
            
def ACEScct_to_linAP1(triplet):
    linAP1 = []
    for channel in triplet:
        if channel <= 0.155251141552511:
            linAP1.append((channel - 0.0729055341958355) / 10.5402377416545)
        elif channel < (log(65504.,2.) + 9.72) / 17.52:
            linAP1.append(2. ** (channel*17.52-9.72))
        else:
            linAP1.append(65504.)
    return linAP1

             
# The next block is based on the Python colour-science package.
# The current version of this script is just a home-made proof-of-concept,
# and not for distribution. But if this is rewritten to become a something
# officially supported by ACES, these transformations must be rewritten 
# from source documents or properly credited/licensed.

PLANCKIAN_TABLE_TUVD = namedtuple('PlanckianTable_Tuvdi', ('Ti', 'ui', 'vi',
                                                           'di'))

CCT_MINIMAL = 1000
CCT_MAXIMAL = 100000
CCT_SAMPLES = 10
CCT_CALCULATION_ITERATIONS = 6

ROBERTSON_ISOTEMPERATURE_LINES_DATA = (
    (0, 0.18006, 0.26352, -0.24341),
    (10, 0.18066, 0.26589, -0.25479),
    (20, 0.18133, 0.26846, -0.26876),
    (30, 0.18208, 0.27119, -0.28539),
    (40, 0.18293, 0.27407, -0.30470),
    (50, 0.18388, 0.27709, -0.32675),
    (60, 0.18494, 0.28021, -0.35156),
    (70, 0.18611, 0.28342, -0.37915),
    (80, 0.18740, 0.28668, -0.40955),
    (90, 0.18880, 0.28997, -0.44278),
    (100, 0.19032, 0.29326, -0.47888),
    (125, 0.19462, 0.30141, -0.58204),
    (150, 0.19962, 0.30921, -0.70471),
    (175, 0.20525, 0.31647, -0.84901),
    (200, 0.21142, 0.32312, -1.0182),
    (225, 0.21807, 0.32909, -1.2168),
    (250, 0.22511, 0.33439, -1.4512),
    (275, 0.23247, 0.33904, -1.7298),
    (300, 0.24010, 0.34308, -2.0637),
    (325, 0.24792, 0.34655, -2.4681),  # 0.24702 ---> 0.24792 Bruce Lindbloom
    (350, 0.25591, 0.34951, -2.9641),
    (375, 0.26400, 0.35200, -3.5814),
    (400, 0.27218, 0.35407, -4.3633),
    (425, 0.28039, 0.35577, -5.3762),
    (450, 0.28863, 0.35714, -6.7262),
    (475, 0.29685, 0.35823, -8.5955),
    (500, 0.30505, 0.35907, -11.324),
    (525, 0.31320, 0.35968, -15.628),
    (550, 0.32129, 0.36011, -23.325),
    (575, 0.32931, 0.36038, -40.770),
    (600, 0.33724, 0.36051, -116.45))
    
ROBERTSON_ISOTEMPERATURE_LINES_RUVT = namedtuple('WyszeckiRobertson_ruvt',
                                                 ('r', 'u', 'v', 't'))

ROBERTSON_ISOTEMPERATURE_LINES = [
    ROBERTSON_ISOTEMPERATURE_LINES_RUVT(*x)
    for x in ROBERTSON_ISOTEMPERATURE_LINES_DATA
]


def uv_to_CCT(uv):
    """
    Returns the correlated colour temperature :math:`T_{cp}` and
    :math:`\Delta_{uv}` from given *CIE UCS* colourspace *uv* chromaticity
    coordinates using *Roberston (1968)* method.

    Parameters
    ----------
    uv : array_like
        *CIE UCS* colourspace *uv* chromaticity coordinates.

    Returns
    -------
    ndarray
        Correlated colour temperature :math:`T_{cp}`, :math:`\Delta_{uv}`.

    References
    ----------
    -   :cite:`AdobeSystems2013`
    -   :cite:`Wyszecki2000y`

    Examples
    --------
    >>> uv = np.array([0.193741375998230, 0.315221043940594])
    >>> uv_to_CCT_Robertson1968(uv)  # doctest: +ELLIPSIS
    array([  6.5000162...e+03,   8.3333289...e-03])
    """

    u, v = uv

    last_dt = last_dv = last_du = 0

    for i in range(1, 31):
        wr_ruvt = ROBERTSON_ISOTEMPERATURE_LINES[i]
        wr_ruvt_previous = ROBERTSON_ISOTEMPERATURE_LINES[i - 1]

        du = 1
        dv = wr_ruvt.t

        length = np.hypot(1, dv)

        du /= length
        dv /= length

        uu = u - wr_ruvt.u
        vv = v - wr_ruvt.v

        dt = -uu * dv + vv * du

        if dt <= 0 or i == 30:
            if dt > 0:
                dt = 0

            dt = -dt

            f = 0 if i == 1 else dt / (last_dt + dt)

            T = 1.0e6 / (wr_ruvt_previous.r * f + wr_ruvt.r * (1 - f))

            uu = u - (wr_ruvt_previous.u * f + wr_ruvt.u * (1 - f))
            vv = v - (wr_ruvt_previous.v * f + wr_ruvt.v * (1 - f))

            du = du * (1 - f) + last_du * f
            dv = dv * (1 - f) + last_dv * f

            length = np.hypot(du, dv)

            du /= length
            dv /= length

            D_uv = uu * du + vv * dv

            break

        last_dt = dt
        last_du = du
        last_dv = dv

    return [T, -D_uv]


def CCT_to_uv(CCT, D_uv=0):
    """
    Returns the *CIE UCS* colourspace *uv* chromaticity coordinates from given
    correlated colour temperature :math:`T_{cp}` and :math:`\Delta_{uv}` using
    *Roberston (1968)* method.

    Parameters
    ----------
    CCT : numeric
        Correlated colour temperature :math:`T_{cp}`.
    D_uv : numeric
        :math:`\Delta_{uv}`.

    Returns
    -------
    ndarray
        *CIE UCS* colourspace *uv* chromaticity coordinates.

    References
    ----------
    -   :cite:`AdobeSystems2013a`
    -   :cite:`Wyszecki2000y`

    Examples
    --------
    >>> CCT = 6500.0081378199056
    >>> D_uv = 0.008333331244225
    >>> CCT_to_uv_Robertson1968(CCT, D_uv)  # doctest: +ELLIPSIS
    array([ 0.1937413...,  0.3152210...])
    """

    r = 1.0e6 / CCT

    for i in range(30):
        wr_ruvt = ROBERTSON_ISOTEMPERATURE_LINES[i]
        wr_ruvt_next = ROBERTSON_ISOTEMPERATURE_LINES[i + 1]

        if r < wr_ruvt_next.r or i == 29:
            f = (wr_ruvt_next.r - r) / (wr_ruvt_next.r - wr_ruvt.r)

            u = wr_ruvt.u * f + wr_ruvt_next.u * (1 - f)
            v = wr_ruvt.v * f + wr_ruvt_next.v * (1 - f)

            uu1 = uu2 = 1.0
            vv1, vv2 = wr_ruvt.t, wr_ruvt_next.t

            length1 = np.hypot(1, vv1)
            length2 = np.hypot(1, vv2)

            uu1 /= length1
            vv1 /= length1

            uu2 /= length2
            vv2 /= length2

            uu3 = uu1 * f + uu2 * (1 - f)
            vv3 = vv1 * f + vv2 * (1 - f)

            len3 = np.sqrt(uu3 * uu3 + vv3 * vv3)

            uu3 /= len3
            vv3 /= len3

            u += uu3 * -D_uv
            v += vv3 * -D_uv

            return [u, v]

#This is the end of the block based on the Python colour-science package.

# This function turns user input into x,y coordinates, 
# even if the user is not careful about syntax.
def userXY_to_xy (user_xy):

    #The user may have included extra spaces.
    xy = str(user_xy).strip()
    xy = re.sub( '\s+', ' ', xy )

    #The user may have omitted commas in lists.       
    xy_string = ''
    for st in range( len(xy) ):
        if st == 0 or st == len(xy) - 1:
            xy_string += xy[st]
        elif xy[st] == ' ' and xy_string[-1].isdigit():
            if xy[st+1].isdigit() or xy[st+1]=='.' or xy[st+1]=='-':
                xy_string += (',')
        else:
            xy_string += xy[st]
                        
    #The user may have omitted parentheses.
    parenths = False
    for st in xy:
        if st == '(' or st == '[':
            parenths == True
            break
    if parenths == False:
        xy = ( '(' + xy_string + ')' )
    else:
        xy = xy_string
    
    #Evaluate the user string, forcing a list even if it evaluates as tuple.     
    return list(eval(xy))

#Convert the user's white balance input to usable CCT and Tint.
def userWB_to_wb(user_cct, user_tint):

    if type(user_cct) == str:
        cct = ''
        for c in user_cct:
            if c.isdigit() or c == '.' or c == '-' : cct += c
        cct = eval(cct)
    else:
        cct = float(user_cct)
        
    if type(user_tint) == str:
        tint = ''
        for c in user_tint:
            if c.isdigit() or c == '.' or c == '-' : tint += c
        tint = eval(tint)
    else:
        tint = float(user_tint)
        
    return [cct,tint]
    


# The main calculator function takes the user's input and returns the targets.
# Even though the chips are semi-hardcoded, they're treated here as user input,
# for future versions in which the user may be able to define the color chips.
def calc_data_set (illuminant, cct, tint, colorChips, colorspace):

    #Convert user entered strings to lists.
    user_illum = userXY_to_xy(illuminant)
    user_whiteBal =userWB_to_wb(cct , tint)
    
    #Get the x,y coordinates for the white balance.
    whiteBal = CCT_to_uv(user_whiteBal[0],user_whiteBal[1])
    whiteBal = uv_to_xy(whiteBal)
        
    #Poynton has given the color chips in Illuminant C, 
    #so we'll need the x,y coordinate for illuminant C.
    IllumC = [0.3101, 0.3163]
    
    #The AP0 x,y white point in given in the ACES documentation:
    AP0 = [0.32168, 0.33767]
    
    #To do chromatic adaptation as gains, we'll need to convert to CAT02.
    
    # We do not add any of what Fairchild calls "discounting," because this has to work 
    # simultaneously for all scene objects (both reflective and emissive), and humans are 
    # most tuned in to nuances of reflective colors. So we'll concentrate on those.
    illum = np.array(xy_to_CAT02(user_illum))
    whiteBal = np.array(xy_to_CAT02(whiteBal))
    IllumC = np.array(xy_to_CAT02(IllumC))
    AP0 = np.array(xy_to_CAT02(AP0))
    
    # Find the white balance ratio, which is the ratio between
    # the camera white balance setting and the actual color of the light.
    wbRatio = illum / whiteBal
    
    #Build the data set output triplets.
    dataSet = []
    for chip in colorChips:
        c = deepcopy(chip[1])
        c = Yxy_to_XYZ(c)
        c = np.array(XYZ_to_CAT02(c))
        
        #Convert reference white from Illuminant C to AP0.
        c *= AP0 / IllumC
        
        # Account for the shift due to the fact that the actual illuminant (as seen through 
        # the lens) may not have been the same color as the user's white balance setting.
        c *= wbRatio
        
        #We're done with the CAT02 adaptation and ready to go to ACES.
        c = CAT02_to_XYZ(c)
        
        c = XYZ_to_linAP0(c)
        
        #Here, we scale the intensity for no reason other than to match the published
        #ACES document.
        c = [x*.202/.198 for x in c]
        
        #Convert to the user's preferred color space.
        if colorspace == 'AP0' :
            dataSet.append([chip[0], c])
        else : 
            c = linAP0_to_linAP1(c)
            if colorspace == 'AP1': dataSet.append([chip[0], c])
            elif colorspace == 'CC': dataSet.append([chip[0], linAP1_to_ACEScc(c)])
            else : dataSet.append([chip[0], linAP1_to_ACEScct(c)])

    return [dataSet, list(user_illum), list(user_whiteBal)]
    
#Tee up special characters for the user interface.
delta = (u'\u0394')
deg = (u'\u00BA')
    
def user_message():
    print ('\n\n-----------------------------------'
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
    
    print ('-----------------------------------\n\n\n'
         'Enter the CIE 1931 x,y coordinates of the illuminant that fell on the Macbeth chart.\n'
         'Mesure the illumant through the same lens you used to shoot the chart, so that\n'
         'any color shift from the lens can be cancelled out in this calculation.\n\n'
         'x,y coordinant of illuminant, through the lens?')
    illuminant = input()

    print ('\nEnter the camera white balance CCT in degrees Kelvin that you are testing\n'
        'or for which you\'re trying to build an IDT.\n'
        'User CCT setting?')
    cct = input()

    print ('\nEnter the camera white balnce tint in ' , end='')
    try: print (delta , end='')
    except: print ('Delta-' , end='')
    print ('uv that you are testing or for which you\'re trying to build an IDT.\n'
        'User tint setting?')
    tint = input()

    print ('\nDo you want the targets for the whole Macbeth chart or just for\n'    
        'middle gray, red, green, blue, cyan, magenta, and yellow?\n'    
        '\t1. 18% gray and Whole Macbeth chart\n'    
        '\t2. 18% gray and Macbeth Gray + RGBCMY')
    chips = input()
    chips = chips.strip()
    
    print ('\nEnter "n" if you would like results as a Nuke node.\n'
        'otherwise, results will be given as human-readable list.')
    
    nuke = input()
    nuke = nuke.strip()
    if nuke.lower() == 'n':
        nuke = True
    else:
        nuke = False
        
    # These are the Macbeth color chips as published by Mr Charles Poynton.
    # I have measured my own personal Macbeth chart and it is incredibly close,
    # so we know those charts are made with impressive reliability and repeatability.
    macbeth =[
    ['18% gray', [.18, 0.31, 0.316]],
    ['dark skin', [0.101, 0.4, 0.35]],
    ['light skin', [0.358, 0.377, 0.345]],
    ['blue sky', [0.193, 0.247, 0.251]],
    ['foliage', [0.133, 0.337, 0.422]],
    ['blue flower', [0.243, 0.265, 0.24]],
    ['bluish green', [0.431, 0.261, 0.343]],
    ['orange', [0.301, 0.506, 0.407]],
    ['purplish blue', [0.12, 0.211, 0.175]],
    ['moderate red', [0.198, 0.453, 0.306]],
    ['purple', [0.066, 0.285, 0.202]],
    ['yellow green', [0.443, 0.38, 0.489]],
    ['orange yellow', [0.431, 0.473, 0.438]],
    ['blue', [0.061, 0.187, 0.129]],
    ['green', [0.234, 0.305, 0.478]],
    ['red', [0.12, 0.539, 0.313]],
    ['yellow', [0.591, 0.448, 0.47]],
    ['magenta', [0.198, 0.364, 0.233]],
    ['cyan', [0.198, 0.196, 0.252]],
    ['white', [0.9, 0.31, 0.316]],
    ['neutral_8', [0.591, 0.31, 0.316]],
    ['neutral_6.5', [0.362, 0.31, 0.316]],
    ['neutral_5', [0.198, 0.31, 0.316]],
    ['neutral_3.5', [0.09, 0.31, 0.316]],
    ['black', [0.031, 0.31, 0.316]],
    ]
    
    #Same data, but from the Python colour-science module instead of from Poynton.
    #Apparently given in D50, per the documentation.
#     COLORCHECKER_2005_DATA = (
#     (1, 'dark skin', np.array([0.4316, 0.3777, 0.1008])),
#     (2, 'light skin', np.array([0.4197, 0.3744, 0.3495])),
#     (3, 'blue sky', np.array([0.2760, 0.3016, 0.1836])),
#     (4, 'foliage', np.array([0.3703, 0.4499, 0.1325])),
#     (5, 'blue flower', np.array([0.2999, 0.2856, 0.2304])),
#     (6, 'bluish green', np.array([0.2848, 0.3911, 0.4178])),
#     (7, 'orange', np.array([0.5295, 0.4055, 0.3118])),
#     (8, 'purplish blue', np.array([0.2305, 0.2106, 0.1126])),
#     (9, 'moderate red', np.array([0.5012, 0.3273, 0.1938])),
#     (10, 'purple', np.array([0.3319, 0.2482, 0.0637])),
#     (11, 'yellow green', np.array([0.3984, 0.5008, 0.4446])),
#     (12, 'orange yellow', np.array([0.4957, 0.4427, 0.4357])),
#     (13, 'blue', np.array([0.2018, 0.1692, 0.0575])),
#     (14, 'green', np.array([0.3253, 0.5032, 0.2318])),
#     (15, 'red', np.array([0.5686, 0.3303, 0.1257])),
#     (16, 'yellow', np.array([0.4697, 0.4734, 0.5981])),
#     (17, 'magenta', np.array([0.4159, 0.2688, 0.2009])),
#     (18, 'cyan', np.array([0.2131, 0.3023, 0.1930])),
#     (19, 'white 9.5 (.05 D)', np.array([0.3469, 0.3608, 0.9131])),
#     (20, 'neutral 8 (.23 D)', np.array([0.3440, 0.3584, 0.5894])),
#     (21, 'neutral 6.5 (.44 D)', np.array([0.3432, 0.3581, 0.3632])),
#     (22, 'neutral 5 (.70 D)', np.array([0.3446, 0.3579, 0.1915])),
#     (23, 'neutral 3.5 (1.05 D)', np.array([0.3401, 0.3548, 0.0883])),
#     (24, 'black 2 (1.5 D)', np.array([0.3406, 0.3537, 0.0311])))

    if chips[0] == '2': colorChips = [macbeth[0]] + [macbeth[22]] + macbeth[13:19]
    else: colorChips = macbeth
    
    print ('\nIn what format would you like your ACES values returned?\n'
        '    1. ACES AP0 linear\n'
        '    2. ACES AP1 linear\n'
        '    3. ACEScc log\n'
        '    4. ACEScct log\n')
    ACES_format = input()
    ACES_format = ACES_format.strip()

    if ACES_format[0] == '2': colorspace = 'AP1'
    elif ACES_format[0] == '3': colorspace = 'CC'
    elif ACES_format[0] == '4': colorspace = 'CCT'
    else: colorspace = 'AP0'

    #Perform the main function.
    dataSet = calc_data_set(illuminant, cct, tint, colorChips, colorspace)

    #Print the results to the screen either as a Nuke Constant node or in a form
    #which is both human-readable and usable as comma-separated.
    print('\n\n-----------------------------------\n'
        'A properly exposed Macbeth chart illuminated with CIE 1931 x,y chromaticity of:\n'
        'x = ' + '{0:.4f}'.format(dataSet[1][0]) + ' , y = ' + '{0:.4f}'.format(dataSet[1][1]) + '\n'
        'as seen through the lens\n'
        'should yield the below RGB triplets if the user White Balance is set to:\n'
        '{0:.0f}'.format(dataSet[2][0]) + 'K , ' + '{0:.4f}'.format(dataSet[2][1]) + ' ' , end='')
    try : print( delta , end='')
    except: print ( 'Delta-' , end='' )
    print('uv.\n\n'
        'Values given in ' , end='')

    if colorspace == 'AP1' : print ('ACES rgb (linear AP1):\n')
    elif colorspace == 'CC' : print ('ACEScc rgb (logarithmic):\n')
    elif colorspace == 'CCT' : print ('ACEScct rgb (logarithmic):\n')
    else : print ('ACES (linear AP0):\n')
    if nuke:
    
        print('\nGiven here as a Nuke Constant node with results as sequential frames:\n')
        
        rCurve = ''
        gCurve = ''
        bCurve = ''
        for chip in dataSet[0]:
            rCurve += (str(chip[1][0]) + ' ')
            gCurve += (str(chip[1][1]) + ' ')
            bCurve += (str(chip[1][2]) + ' ')
        
        nukeNode = 'Constant {\n inputs 0\n channels rgb\n color {{curve '
        nukeNode += rCurve
        nukeNode += '} {curve '
        nukeNode += gCurve
        nukeNode += '} {curve '
        nukeNode += bCurve
        nukeNode += '} {curve 0}}\n format "1920 1080 0 0 1920 1080 1 HD_1080"\n name ACES_targets\n selected true\n}\n\n'
        
        print (nukeNode)
            
        
    else:  
        for chip in dataSet[0]:
            print (chip[0] + ' , ' , end = '')
            for c in range(3):
                if c<2: print ('{0:.3f}'.format(chip[1][c]) , end=' , ')
                else: print ('{0:.3f}'.format(chip[1][c]) , end='\n' )
        print ('-----------------------------------\n')
    
    

user_message()

#loop the interface
looper = True
while looper == True:
    try:
        user_interface()
        goAgain = input('Enter "q" to quit or any other key to run again.')
    except:
        goAgain = input('\n\nError: could not interpret input. "q" to quit or any other key to try again.')
    if goAgain.lower() == 'q' or goAgain.lower() == 'quit':
        looper = False


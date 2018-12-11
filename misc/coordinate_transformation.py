#
# -*- coding: iso-8859-1 -*-
###########################################################################
# 
# File:            coordinates.py
#
# Author:          Olli Lammi
#
# Version:         0.7a
#
# Date:            25.02.2008
#
# Functions:       KKJxy_to_WGS84lalo
#                  WGS84lalo_to_KKJxy
#                  KKJxy_to_KKJlalo
#                  KKJlalo_to_KKJxy
#                  KKJlalo_to_WGS84lalo
#                  WGS84lalo_to_KKJlalo
#                  KKJ_Zone_I
#                  KKJ_Zone_Lo
#                  WGS84lalo_to_GoogleMapsXY
#                  Str_to_CoordinateValue
#                  KKJxy_in_Finland
#                   
# Description:     Coordinate system functions. Initial code and alhorithms
#                  extracted from Viestikallio's PHP-source code 
#                  (http://www.viestikallio.fi/tools/kkj-wgs84.php). 
#                  Python translation and changes by Olli Lammi.
#
#                  Google Maps functions developed based on knowledge on
#                  multiple sites in the Internet, simplifying the algorithms
#                  and traditional trial error method.
# 
#                  NOTE!: The coordinate functions are developed to work only
#                  with coordinates that are in the area of Finland.
#
# Version history: ** 05.09.2005 v0.1a (Olli Lammi) **
#                  First version. Translated partially from PHP to Python.
#                  Original PHP code source: www.viestikallio.fi
# 
#                  ** 26.09.2005 v0.2a (Olli Lammi) **
#                  Included the WGS84_to_KKJxy -conversion.
#                  Changed all interfaces to work with degrees when
#                  using angle values
#                  Altered the KKJ zone info function to a
#                  lookup dictionary (KKJ_ZONE_INFO).
#                  Added function to calculate propable KKJ band
#                  from KKJ longitude (KKJ_Zone_Lo).
#                  
#                  ** 01.10.2005 v0.3a (Olli Lammi) **
#                  Small changes to function interfaces.
#
#                  ** 14.10.2005 v0.4a (Olli Lammi) **
#                  Added support for KKJ-bands 0 and 5.
#
#                  ** 26.01.2007 v0.5a (Olli Lammi) **
#                  Added prototype function to convert WGS84-coordinates
#                  to Google Maps URL tile parameter x and y.
#
#                  ** 19.02.2008 v0.6a (Olli Lammi) **
#                  Added funcion to convert WGS84 coordinate string to
#                  scalar coordinate value.
#
#                  ** 25.02.2008 v0.7a (Olli Lammi) **
#                  Added funcion to tell approximately whether given KKJ
#                  koordinate is in the are of Finland.
#
###########################################################################

# Imports

import sys, os, string
import math
import re


###########################################################################

# Constants

# Longitude0 and Center meridian of KKJ bands
KKJ_ZONE_INFO = { 0: (18.0,  500000.0), \
                  1: (21.0, 1500000.0), \
                  2: (24.0, 2500000.0), \
                  3: (27.0, 3500000.0), \
                  4: (30.0, 4500000.0), \
                  5: (33.0, 5500000.0), \
                }


# Functions

###########################################################################
# Function:  KKJxy_to_WGS84lalo
###########################################################################
# Input:     dictionary with ['P'] is KKJ Northing
#                            ['I'] in KKJ Eeasting
# Output:    dictionary with ['La'] is latitude in degrees (WGS84)
#                            ['Lo'] is longitude in degrees (WGS84)
###########################################################################

def KKJxy_to_WGS84lalo(KKJin):  
  KKJz = KKJxy_to_KKJlalo(KKJin)
  WGS = KKJlalo_to_WGS84lalo(KKJz)

  return WGS



###########################################################################
# Function:  WGS84lalo_to_KKJxy
###########################################################################
# Input:     dictionary with ['La'] is latitude in degrees (WGS84)
#                            ['Lo'] is longitude in degrees (WGS84)
# Output:    dictionary with ['P'] is KKJ Northing
#                            ['I'] in KKJ Eeasting
###########################################################################

def WGS84lalo_to_KKJxy(WGSin):
  KKJlalo = WGS84lalo_to_KKJlalo(WGSin);

  ZoneNumber = KKJ_Zone_Lo(KKJlalo['Lo'])
  KKJxy = KKJlalo_to_KKJxy(KKJlalo, ZoneNumber)

  return KKJxy



###########################################################################
# Function:  KKJlalo_to_WGS84lalo
###########################################################################

def KKJlalo_to_WGS84lalo(KKJ):
  La = KKJ['La']
  Lo = KKJ['Lo']

  dLa = math.radians(  0.124867E+01       + \
                      -0.269982E+00 * La + \
                     0.191330E+00 * Lo + \
                     0.356119E-02 * La * La + \
                     -0.122312E-02 * La * Lo + \
                     -0.335514E-03 * Lo * Lo ) / 3600.0
  dLo = math.radians( -0.286111E+02       + \
                    0.114183E+01 * La + \
                    -0.581428E+00 * Lo + \
                    -0.152421E-01 * La * La + \
                    0.118177E-01 * La * Lo + \
                    0.826646E-03 * Lo * Lo ) / 3600.0

  WGS = {}
  WGS['La'] = math.degrees(math.radians(KKJ['La']) + dLa)
  WGS['Lo'] = math.degrees(math.radians(KKJ['Lo']) + dLo)

  return WGS



###########################################################################
# Function:  WGS84lalo_to_KKJlalo
###########################################################################

def WGS84lalo_to_KKJlalo(WGS):
  La = WGS['La']
  Lo = WGS['Lo']

  dLa = math.radians( -0.124766E+01       + \
                    0.269941E+00 * La + \
                    -0.191342E+00 * Lo + \
                    -0.356086E-02 * La * La + \
                    0.122353E-02 * La * Lo + \
                    0.335456E-03 * Lo * Lo ) / 3600.0
  
  dLo = math.radians(  0.286008E+02       + \
                     -0.114139E+01 * La + \
                     0.581329E+00 * Lo + \
                     0.152376E-01 * La * La + \
                     -0.118166E-01 * La * Lo + \
                     -0.826201E-03 * Lo * Lo ) / 3600.0

  KKJ = {}
  KKJ['La'] = math.degrees(math.radians(WGS['La']) + dLa)
  KKJ['Lo'] = math.degrees(math.radians(WGS['Lo']) + dLo)

  return KKJ



###########################################################################
# Function:  KKJxy_to_KKJlalo
###########################################################################

def KKJxy_to_KKJlalo(KKJ):  
  #
  # Scan iteratively the target area, until find matching
  # KKJ coordinate value.  Area is defined with Hayford Ellipsoid.
  #  
  LALO = {}

  ZoneNumber = KKJ_Zone_I(KKJ['I'])
    
  MinLa = math.radians(59.0)
  MaxLa = math.radians(70.5)
  MinLo = math.radians(18.5)
  MaxLo = math.radians(32.0)

  i = 1
  while (i < 35):
    DeltaLa = MaxLa - MinLa
    DeltaLo = MaxLo - MinLo

    LALO['La'] = math.degrees(MinLa + 0.5 * DeltaLa)
    LALO['Lo'] = math.degrees(MinLo + 0.5 * DeltaLo)

    KKJt = KKJlalo_to_KKJxy(LALO, ZoneNumber)

    if (KKJt['P'] < KKJ['P']):
      MinLa = MinLa + 0.45 * DeltaLa
    else:
      MaxLa = MinLa + 0.55 * DeltaLa

    if (KKJt['I'] < KKJ['I']):
      MinLo = MinLo + 0.45 * DeltaLo
    else:
      MaxLo = MinLo + 0.55 * DeltaLo

    i = i + 1

  return LALO



###########################################################################
# Function:  KKJlalo_to_KKJxy
###########################################################################

def KKJlalo_to_KKJxy(INP, ZoneNumber):
  Lo = math.radians(INP['Lo']) - math.radians(KKJ_ZONE_INFO[ZoneNumber][0])

  a  = 6378388.0            # Hayford ellipsoid
  f  = 1/297.0

  b  = (1.0 - f) * a
  bb = b * b              
  c  = (a / b) * a        
  ee = (a * a - bb) / bb  
  n = (a - b)/(a + b)     
  nn = n * n              

  cosLa = math.cos(math.radians(INP['La']))

  NN = ee * cosLa * cosLa 

  LaF = math.atan(math.tan(math.radians(INP['La'])) / math.cos(Lo * math.sqrt(1 + NN)))

  cosLaF = math.cos(LaF)

  t   = (math.tan(Lo) * cosLaF) / math.sqrt(1 + ee * cosLaF * cosLaF)

  A   = a / ( 1 + n )

  A1  = A * (1 + nn / 4 + nn * nn / 64)

  A2  = A * 1.5 * n * (1 - nn / 8)

  A3  = A * 0.9375 * nn * (1 - nn / 4)

  A4  = A * 35/48.0 * nn * n

  OUT = {}
  OUT['P'] = A1 * LaF - \
        A2 * math.sin(2 * LaF) + \
            A3 * math.sin(4 * LaF) - \
                A4 * math.sin(6 * LaF)
  OUT['I'] = c * math.log(t + math.sqrt(1+t*t)) + \
        500000.0 + ZoneNumber * 1000000.0

  return OUT



###########################################################################
# Function:  KKJ_Zone_I
###########################################################################

def KKJ_Zone_I(KKJI):
  ZoneNumber = math.floor((KKJI/1000000.0))
  if ZoneNumber < 0 or ZoneNumber > 5:
      ZoneNumber = -1
      
  return ZoneNumber



###########################################################################
# Function:  KKJ_Zone_Lo
###########################################################################

def KKJ_Zone_Lo(KKJlo):
  # determine the zonenumber from KKJ easting
  # takes KKJ zone which has center meridian
  # longitude nearest (in math value) to
  # the given KKJ longitude
  ZoneNumber = 5
  while ZoneNumber >= 0:
    if math.fabs(KKJlo - KKJ_ZONE_INFO[ZoneNumber][0]) <= 1.5:
      break
    ZoneNumber = ZoneNumber - 1
            
  return ZoneNumber



###########################################################################
# Function:  WGS84lalo_to_GoogleMapsXY
###########################################################################
# Input:     dictionary with ['La'] is latitude in degrees (WGS84)
#                            ['Lo'] is longitude in degrees (WGS84)
#            google zoom factor (integer between 0 and 17, 0 for no zoom 
#                            and 17 for maximum zoom)
# Output:    dictionary with ['x'] is Google maps URL x parameter (tile number)
#                            ['y'] in Google maps URL y parameter (tile number)
###########################################################################

def WGS84lalo_to_GoogleMapsXY(WGSin, zoom):
  # Google maps maximum zoom factor (min = 0)
  MAXZOOM = 17.0

  worldwidth = math.pow(2.0, (MAXZOOM - zoom))
  x = WGSin['Lo'] + 180.0
  y = math.log(math.tan((math.pi / 4.0) + ((0.5 * math.pi * WGSin['La']) / 180.0))) / math.pi
  if y <-0.9999:
    y = -0.9999
  if y > 0.9999:
    y = 0.9999 
  y = (-90.0 * y + 90.0)
  
  out = {}
  out['x'] = (int) (math.floor((x / 360.0) * worldwidth))
  out['y'] = (int) (math.floor((y / 180.0) * worldwidth))

  return out



###########################################################################
# Function:  Str_to_CoordinateValue
###########################################################################
# Input:     string with a coordinate value in some WGS84-format
# 
# Output:    floating point value representing the given coordinate
#            Value returned is INVALID_COORDINATE if the given string cannot be
#            interpreted as WGS84 coordinate value. 
###########################################################################
INVALID_COORDINATE = -99999.99
def Str_to_CoordinateValue(WGSstr):
  
  # case 1: 61,27,4.96 (for 61 degrees, 27 minutes, 4.96 seconds)
  regexp1 = '^(?P<sig>-?)(?P<deg>\d+),(?P<min>\d+),(?P<sec>\d+(\.\d+)?)$'
  mo = re.match(regexp1, WGSstr)
  if (mo != None):
    value = string.atof( mo.group('deg') ) + string.atof( mo.group('min') ) / 60.0 + string.atof( mo.group('sec') ) / 3600.0
    if ( mo.group('sig') == '-' ):
      value = -value
    return value

  # case 2: 61,27.083 (for 61 degrees, 27.083 minutes)
  regexp2 = '^(?P<sig>-?)(?P<deg>\d+),(?P<min>\d+(\.\d+))?$'
  mo = re.match(regexp2, WGSstr)
  if (mo != None):
    value = string.atof( mo.group('deg') ) + string.atof( mo.group('min') ) / 60.0
    if ( mo.group('sig') == '-' ):
      value = -value
    return value

  # case 3: 61.451378 (for 61.451378 degrees)
  regexp3 = '^(?P<sig>-?)(?P<deg>\d+\.\d+)$'
  mo = re.match(regexp3, WGSstr)
  if (mo != None):
    value = string.atof( mo.group('deg') )
    if ( mo.group('sig') == '-' ):
        value = -value
    return value

  # case: other
  return INVALID_COORDINATE 



###########################################################################
# Function:  KKJxy_in_Finland
###########################################################################
# Input:     dictionary with ['P'] is KKJ Northing
#                            ['I'] in KKJ Eeasting
# Output:    truth value telling whether the given coordinate is 
#            _*_*_approximately_*_*_ in Finnish area. 
###########################################################################

FINLAND_AREA_ZONE3_BOXES = [ [6641806, 3069283, 6737806, 3161283], \
                             [6641806, 3161283, 7067806, 3755283], \
                             [7067806, 3253283, 7581806, 3691283], \
                             [7581806, 3413283, 7786097, 3595283], \
                             [7581806, 3219283, 7726097, 3413283]
                            ]

def KKJxy_in_Finland(KKJ):
  try:
    # Move the coordinates to zone 3
    lalo = KKJxy_to_KKJlalo(KKJ)
    xy = KKJlalo_to_KKJxy(lalo, 3)
 
    P = int(round(xy['P']))
    I = int(round(xy['I']))

    # go through the boxes and check if coordinates are inside one of them 
    for box in FINLAND_AREA_ZONE3_BOXES:
      if P >= box[0] and P <= box[2] and I >= box[1] and I <= box[3]: 
        return 1
  except:
    # If KKJ-conversion functions fail, assume it is not valid Finnish location
    return 0
     
  return 0

     
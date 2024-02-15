import math
import numpy as np

def metersToMiles(lengthInMeters):
    return lengthInMeters / 1852.0

def milesToMeters(lengthInMiles):
    return lengthInMiles * 1852.0

def feetToMeters(feet):
    return feet * 0.3048

def degToRad(deg):
    return deg / 180.0 * math.pi

def radToDeg(rad):
    return rad / math.pi * 180.0

earthRadius = 6371000 # meters

# lat/lon/height (deg, m amsl) to ecef (m)
def llh2ecef(llh):
    r = llh[2] + earthRadius
    x = r * math.cos(degToRad(llh[0])) * math.cos(degToRad(llh[1]))
    y = r * math.cos(degToRad(llh[0])) * math.sin(degToRad(llh[1]))
    z = r * math.sin(degToRad(llh[0]))
    return np.array([x, y, z])

# ecef (m) to lat/lon/height (deg, m amsl)
def ecef2llh(ecef):
    if ecef[0] == 0 and ecef[1] == 0:
        if ecef[2] > 0: lat = 90
        elif ecef[2] < 0: lat = -90
        else: lat = 0
    else:
        lat = radToDeg(math.atan(ecef[2] / math.sqrt(ecef[0]*ecef[0]+ecef[1]*ecef[1])))
    lon = radToDeg(math.atan2(ecef[1], ecef[0]))
    height = math.sqrt(ecef[0]*ecef[0]+ecef[1]*ecef[1]+ecef[2]*ecef[2]) - earthRadius
    return np.array([lat, lon, height])
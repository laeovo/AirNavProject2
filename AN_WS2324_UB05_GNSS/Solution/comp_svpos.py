"""
    Exercise 05 - comp_svpos.py
    
    - Calculate the satellite position at a particular time (given)
    - Convert teh position to ENU at a specific origin (given)
    
    Copyright (c) - Maarten Uijt de Haag
"""

import time
import sys
import math
import numpy as np
import matplotlib.pyplot as plt

from ephcal import *
from llh2ecef import *
from ecef2enu import *

# Load the Ephmereris data
ephem = np.loadtxt('ephem.dat')

# Check out what satellites are described in this data set
svid = np.nonzero(ephem[:,0])[0]

# Time at which the signal was transmitted
tot = 307741.994

# Give the origin's coordinates
orgllh = np.array([[52.8085525646001*np.pi/180.0],[13.7592384849522*np.pi/180.0],[0]])
    
# Convert origin to ECEF
orgece = llh2ecef(orgllh)
    
    
# Go through all satellites in Ephemeris data set
el = np.zeros((len(svid),1))
az = np.zeros((len(svid),1))
for ii in range(0,len(svid)) :
    
    # Calculate the satellite position (in ECEF) 
    (svpos, svclock, ecode) = ephcal(tot, ephem, svid[ii])
    
    
    # Convert satellite position to ENU
    svenu = ecef2enu ( svpos, orgece, orgllh) 
    
    # Compute the satellite Elevation
    el[ii] = math.atan2(svenu[2],math.sqrt(svenu[0]**2 + svenu[1]**2))*180.0/np.pi
    az[ii] = math.atan2(svenu[0],svenu[1])*180.0/np.pi
    
    # Print the azimuth and elevation of the satellite
    print('Azimuth: %5.2lf, elevation: %5.2lf' % (az[ii], el[ii]))


"""
    Exercise 05 - comp_tot_svpos_template.py
    
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

# GPS value for speed of light
vlight   = 299792458.0       

# WGS-84 value for earth rotation rate
eanvel = 7.2921151467e-05    

# Load the raw data
data = np.loadtxt('pr.dat')

# Extract satellite ids and typecast to integer
svid = data[:,0].astype(int)

# Extract pseudoranges
pr = data[:,1]
 
# Load the Ephmereris data
ephem = np.loadtxt('ephem.dat')

# Time at which the signal was RECEIVED
tor = 307741.994

# Give the origin's coordinates
orgllh = np.array([[52.8085525646001*np.pi/180.0],[13.7592384849522*np.pi/180.0],[0]])
    
# Convert origin to ECEF
orgecef = llh2ecef(orgllh)

# Go through all satellites in Ephemeris data set
el = np.zeros((len(svid),1))
az = np.zeros((len(svid),1))
for ii in range(0,len(svid)) :
    
    # Implement the time compensation method from the slides (Part 2)
    
    
    # Correct satellite position for earth rotation during
    # travel time from satellite to user GPS receiver (Part 3)
    
    
    # Convert satellite position to ENU
    
    
    # Compute the satellite elevation and azimuth (from Part 1)
  
    
    # Print the azimuth and elevation of the satellite
    print('Azimuth: %5.2lf, elevation: %5.2lf' % (az[ii], el[ii]))


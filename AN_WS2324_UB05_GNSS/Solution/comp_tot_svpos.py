"""
    Exercise 05 - comp_tot_svpos.py
    
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
orgece = llh2ecef(orgllh)

# Go through all satellites in Ephemeris data set
el = np.zeros((len(svid),1))
az = np.zeros((len(svid),1))
for ii in range(0,len(svid)) :
    
    # Implement the time compensation method from the slides
    tot = tor - pr[ii]/vlight
    
    # Calculate the satellite position (in ECEF) 
    (svpos, svclock, ecode) = ephcal(tot, ephem, svid[ii])
    
    if ( abs(svclock) < 0.01) :
        
        tot = tor - pr[ii]/vlight - svclock
        
        (svpos, svclock, ecode) = ephcal(tot, ephem, svid[ii])
    
    # Correct satellite position for earth rotation during
    # travel time from satellite to user GPS receiver
    # Note: it is not necessary to correct for the receiver
    # clock offset since this particular GPS receiver
    # synchronizes its measurements to GPS time such that the
    # receiver clock offset is negligible (< 100 ns)
    travel_time = pr[ii] / vlight;
    if ( abs(svclock) < 0.01) :
        travel_time = travel_time + svclock;
    else :
        print('Crap')
    
    # Make sure that the travel_time is not unreasonable
    if (travel_time < 0.06) :
        travel_time = 0.06
        print('Crap 2')
    if (travel_time > 0.10) :
        travel_time = 0.10
        print('Crap 3')
        
    alfa = travel_time * eanvel;
    svecef = 0.0*svpos
    svecef[0] = svpos[0] + alfa * svpos[1];
    svecef[1] = svpos[1] - alfa * svpos[0];
    svecef[2] = svpos[2]
    
   
    
    # Convert satellite position to ENU
    svenu = ecef2enu ( svecef, orgece, orgllh) 
    
    # Compute the satellite elevation and azimuth
    el[ii] = math.atan2(svenu[2],math.sqrt(svenu[0]**2 + svenu[1]**2))*180.0/np.pi
    az[ii] = math.atan2(svenu[0],svenu[1])*180.0/np.pi
    
    # Print the azimuth and elevation of the satellite
    print('Azimuth: %5.2lf, elevation: %5.2lf' % (az[ii], el[ii]))


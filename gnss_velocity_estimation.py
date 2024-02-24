import numpy as np
from ephcal import *
import math
from matplotlib import pyplot as plt
from llh2ecef import *

# GPS value for speed of light
vlight = 299792458.0 # m/s

# WGS-84 value for earth rotation rate
eanvel = 7.2921151467e-05 # rad/s

raw = np.loadtxt("group01_raw.dat")
raw = raw[raw[:, 1] == 0] # only GPS data
raw = raw[raw[:, 5] == 0] # only L1 band
raw = raw[raw[:, 8] >= 3] # only data points with tracking status >= 3

# make sure the measurements are sorted by time
for i in range(len(raw)-1):
    if raw[i+1, 0] - raw[i, 0] >= 0: continue
    print("Measurements are not in timely order")
    exit(1)

ephem = np.loadtxt("group01_ephem.dat")

previous_tor = -1
current_tor = 0
satellites_previous = []
satellites_current = []
for raw_row in range(len(raw)):
    if raw[raw_row, 0] == current_tor: continue
    previous_tor = current_tor
    current_tor = raw_row[raw_row, 0]

    # identify satellites common between previous and current tor
    satellites_previous = satellites_current
    satellites_current = raw[raw[:, 0] == current_tor][:,2]
    common_satellite_ids = np.intersect1d(satellites_previous, satellites_current)

    # compute single difference between common satellites
    delta_Phi = np.zeros(len(common_satellite_ids))
    for i in range(len(common_satellite_ids)):
        delta_Phi[i] = raw[raw[:, 0] == current_tor and raw[:,2] == common_satellite_ids[i]] - raw[raw[:, 0] == previous_tor and raw[:,2] == common_satellite_ids[i]]

    # calculate unit vectors to satellites at current time
    e_current = np.zeros((len(common_satellite_ids), 3))
    for i in range(len(common_satellite_ids)):
        e_current[i,:] = 

plt.scatter(positions[:,2], positions[:,1], c=positions[:,3])
plt.colorbar()
plt.xlabel("Longitude (degrees)")
plt.ylabel("Latitude (degrees)")
plt.title("Estimated positions in lat/lon/height (deg/deg/m)")
plt.show()


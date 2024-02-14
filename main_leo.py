import numpy as np
from ephcal import *

# GPS value for speed of light
vlight   = 299792458.0

raw = np.loadtxt("group01_raw.dat")
raw = raw[raw[:, 1] == 0] # only GPS data
raw = raw[raw[:, 5] == 0] # only L1 band
raw = raw[raw[:, 8] >= 3] # only data points with tracking status >= 3

ephem = np.loadtxt("group01_ephem.dat")



for tor in np.unique(raw[:,0]):
    data_points = raw[raw[:, 0] == tor]
    H = np.array([])
    for row in range(len(data_points)):
        datarow = data_points[row]
        current_gps_time = datarow[0]
        svid = datarow[2]
        pr = datarow[3]

        # compute time of transmission (tot)
        tot_first_guess = current_gps_time - pr / vlight
        (svpos, svclock, ecode) = ephcal(tot_first_guess, ephem, svid)
        updated_tot = tot_first_guess - svclock

        # with new tot, compute position of satellite
        (svpos, svclock, ecode) = ephcal(updated_tot, ephem, svid)
        svpos = np.concatenate((1/np.linalg.norm(svpos)*svpos, np.array([[1]])), axis=0).T

        if len(H) == 0: H = np.array(svpos)
        else: H = np.concatenate((H, svpos), axis=0)

    if len(H) >= 4:
        # now we can compute an actual position, because we see at least four satellites
        



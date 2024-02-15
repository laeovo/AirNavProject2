import numpy as np
from ephcal import *
import math

# GPS value for speed of light
vlight = 299792458.0

raw = np.loadtxt("group01_raw.dat")
raw = raw[raw[:, 1] == 0] # only GPS data
raw = raw[raw[:, 5] == 0] # only L1 band
raw = raw[raw[:, 8] >= 3] # only data points with tracking status >= 3

ephem = np.loadtxt("group01_ephem.dat")

positions = np.array([])

print(len(np.unique(raw[:,0])), "different times")

for tor in np.unique(raw[:,0]):
    print("Time:", tor)
    data_points = raw[raw[:, 0] == tor]
    satellites = np.array([])
    z = np.array([])
    for row in range(len(data_points)):
        datarow = data_points[row]
        current_gps_time = tor
        svid = datarow[2]
        pr = datarow[3]

        # compute time of transmission (tot)
        tot_first_guess = current_gps_time - pr / vlight
        (svpos, svclock, ecode) = ephcal(tot_first_guess, ephem, svid)
        tot_updated = tot_first_guess - svclock

        # with new tot, compute position of satellite
        (svpos, svclock, ecode) = ephcal(tot_updated, ephem, svid)

        if len(satellites) == 0:
            satellites = np.array(svpos.T)
            z = np.array([pr])
        else:
            satellites = np.concatenate((satellites, svpos.T), axis=0)
            z = np.concatenate((z, [pr]), axis=0)

    if len(satellites) < 4: continue # need at least four satellites to determine user position

    # now that we see at least four satellites, we can compute an actual user position
    x_hat = np.array([0, 0, 0])

    size_of_last_step = 100000 # arbitrary large value to get the while loop going
    counter = 0
    while size_of_last_step > 0.3 and counter < 2000:
        counter += 1
        # print("      x_hat =", x_hat, ", distance to center =", np.linalg.norm(x_hat))
        # what would we measure if we were at point x_hat and didn't have any errors?
        z_hat = np.zeros(satellites.shape[0])
        for i in range(satellites.shape[0]):
            z_hat[i] = np.linalg.norm(x_hat - satellites[i])

        # print("We actually measured:           ", z[0:3], "...")
        # print("At x_hat we would have measured:", z_hat[0:3], "...")

        # how much are we off from what we actually measured?
        delta_z_hat = z - z_hat

        # compute H matrix based on current estimate
        H = np.zeros(satellites.shape)
        for i in range(len(H)): H[i] = (x_hat - satellites[i]) / np.linalg.norm(x_hat - satellites[i])

        # compute step into the "right" direction
        delta_x_hat = 0.5 * np.linalg.inv(H.T @ H) @ H.T @ delta_z_hat
        size_of_last_step = np.linalg.norm(delta_x_hat)
        # print("delta_x_hat =", delta_x_hat, ", length =", size_of_last_step)

        # update estimate
        x_hat = x_hat + delta_x_hat

    positions = np.concatenate((positions, [tor, x_hat]), axis=0)



print(positions)


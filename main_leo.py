import numpy as np
from ephcal import *
import math
from matplotlib import pyplot as plt
from llh2ecef import *

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
    x_hat = np.array([0, 0, 0, 0])

    size_of_last_step = 100000 # arbitrary large value to get the while loop going
    counter = 0
    while size_of_last_step > 0.3 and counter < 2000:
        counter += 1
        # print("      x_hat =", x_hat, ", distance to center =", np.linalg.norm(x_hat))
        # what would we measure if we were at point x_hat and didn't have any errors?
        z_hat = np.zeros(satellites.shape[0])
        for i in range(satellites.shape[0]):
            z_hat[i] = np.linalg.norm(x_hat[0:3] - satellites[i])

        # print("We actually measured:           ", z[0:3], "...")
        # print("At x_hat we would have measured:", z_hat[0:3], "...")

        # how much are we off from what we actually measured?
        delta_z_hat = z - z_hat

        # compute H matrix based on current estimate
        H = np.zeros((satellites.shape[0], 4))
        for i in range(len(H)): H[i] = np.concatenate(((x_hat[0:3] - satellites[i]) / np.linalg.norm(x_hat[0:3] - satellites[i]), [1]), axis=0)

        # compute step into the "right" direction
        delta_x_hat = np.linalg.inv(H.T @ H) @ H.T @ delta_z_hat
        size_of_last_step = np.linalg.norm(delta_x_hat[0:3])
        # print("delta_x_hat =", delta_x_hat, ", length =", size_of_last_step)

        # update estimate
        x_hat = x_hat + delta_x_hat

    new_position = np.zeros((1, 4))
    new_position[0, 0] = tor
    new_position[0, 1:4] = ecef2llh(x_hat[0:3])
    if new_position[0, 3] < -90000: continue
    if len(positions) == 0: positions = new_position
    else: positions = np.concatenate((positions, new_position), axis=0)
    # positions = np.concatenate((positions, [np.concatenate(([tor], [x_hat[0:3]]), axis=1)]), axis=0)


print(positions)
# plt.scatter(positions[:,0], positions[:,3], c=positions[:,3])
plt.plot(positions[:,0], positions[:,3])
plt.show()


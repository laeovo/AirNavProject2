import numpy as np
from ephcal import *
import math

# GPS value for speed of light
vlight = 299792458.0 # m/s

# WGS-84 value for earth rotation rate
eanvel = 7.2921151467e-05 # rad/s

# calculates space vehicle position, given time-of-reception, single pseudorange, SV ID and ephem data
def calculate_sv_pos(time_of_reception, pr, svid, ephem):
    # compute time of transmission (tot)
    tot_first_guess = time_of_reception - pr / vlight
    (svpos, svclock, ecode) = ephcal(tot_first_guess, ephem, svid)
    tot_updated = tot_first_guess - svclock

    # with new tot, compute position of satellite
    (svpos, svclock, ecode) = ephcal(tot_updated, ephem, svid)

    # apply svclock error to pseudorange
    pr += svclock * vlight

    # correct for Sagnac effect
    delta_t = pr / vlight
    alpha = delta_t * eanvel
    svpos = np.array([[math.cos(alpha), math.sin(alpha), 0],
                      [-math.sin(alpha), math.cos(alpha), 0],
                      [0, 0, 1]]) @ svpos

    return svpos.T, pr

# calculates user position in ECEF and user clock error, given pseudoranges and sv positions
def calculate_least_squares_solution(pseudoranges, satellitePositions):
    # x_hat is the current position estimate; last component is the distance error produced by the user clock (c * delta t_u)
    x_hat = np.array([0, 0, 0, 0])

    size_of_last_step = 100000  # arbitrary large value to get the while loop going
    counter = 0
    while size_of_last_step > 0.3 and counter < 2000:
        counter += 1

        z_hat = np.zeros(satellitePositions.shape[0])
        for i in range(satellitePositions.shape[0]):
            z_hat[i] = np.linalg.norm(x_hat[0:3] - satellitePositions[i]) + x_hat[3]

        # how much are we off from what we actually measured?
        delta_z_hat = pseudoranges - z_hat

        # compute H matrix based on current estimate
        H = np.zeros((satellitePositions.shape[0], 4))
        for i in range(len(H)):
            H[i, 0:3] = (x_hat[0:3] - satellitePositions[i]) / np.linalg.norm(x_hat[0:3] - satellitePositions[i])
            H[i, 3] = 1

        # compute step into the "right" direction
        delta_x_hat = np.linalg.inv(H.T @ H) @ H.T @ delta_z_hat
        size_of_last_step = np.linalg.norm(delta_x_hat)

        # update estimate
        x_hat = x_hat + delta_x_hat

    return x_hat[0:3], x_hat[3]
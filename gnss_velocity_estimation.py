import numpy as np
from ephcal import *
import math
from matplotlib import pyplot as plt
from llh2ecef import *
from gnss_position_estimation import calculate_sv_pos, calculate_least_squares_solution

# GPS value for speed of light
vlight = 299792458.0 # m/s

# WGS-84 value for earth rotation rate
eanvel = 7.2921151467e-05 # rad/s

raw = np.loadtxt("group01_raw.dat")
raw = raw[raw[:, 1] == 0] # only GPS data
raw = raw[raw[:, 5] == 0] # only L1 band
raw = raw[raw[:, 8] >= 3] # only data points with tracking status >= 3

# make sure the measurements are sorted by time
for j in range(len(raw) - 1):
    if raw[j + 1, 0] - raw[j, 0] >= 0: continue
    print("Measurements are not in timely order")
    exit(1)

ephem_from_file = np.loadtxt("group01_ephem.dat")

previous_tor = -1
current_tor = 0
satellite_ids_previous = []
satellite_ids_current = []
velocities = np.zeros((len(np.unique(raw[:, 0]))-1, 3))
times = np.zeros(len(np.unique(raw[:, 0]))-1)
time_counter = 0
for raw_row in range(len(raw)):
    if raw[raw_row, 0] == current_tor: continue
    previous_tor = current_tor
    current_tor = raw[raw_row, 0]
    measurements_previous = raw[raw[:, 0] == previous_tor]
    measurements_current = raw[raw[:, 0] == current_tor]

    # step 1: identify satellites common between previous and current tor
    satellite_ids_previous = satellite_ids_current
    satellite_ids_current = raw[raw[:, 0] == current_tor][:, 2]
    common_satellite_ids = np.intersect1d(satellite_ids_previous, satellite_ids_current)
    if len(common_satellite_ids) < 4: continue # not enough common satellites

    # step 2: compute single difference between common satellites
    delta_Phi = np.zeros(len(common_satellite_ids))
    for j in range(len(common_satellite_ids)):
        svid = common_satellite_ids[j]
        Phi_previous = measurements_previous[measurements_previous[:, 2] == svid][0, 4]
        Phi_current = measurements_current[measurements_current[:, 2] == svid][0, 4]
        delta_Phi[j] = Phi_current - Phi_previous

    # step 3: calculate unit vectors to satellites at current time
    satellite_positions_current = np.zeros((len(common_satellite_ids), 3))
    updated_pseudoranges_current = np.zeros(len(common_satellite_ids))
    for j in range(len(common_satellite_ids)):
        svid = common_satellite_ids[j]
        old_pseudorange = measurements_current[measurements_current[:, 2] == svid][0, 3]
        sv_position, updated_pseudorange = calculate_sv_pos(current_tor, old_pseudorange, svid, ephem_from_file)
        satellite_positions_current[j] = sv_position
        updated_pseudoranges_current[j] = updated_pseudorange
    user_position_current_ecef, user_clock_error = calculate_least_squares_solution(updated_pseudoranges_current, satellite_positions_current)
    e_current = np.zeros((len(common_satellite_ids), 3))
    for j in range(len(common_satellite_ids)):
        e_current[j, :] = (satellite_positions_current[j] - user_position_current_ecef) / np.linalg.norm(satellite_positions_current[j] - user_position_current_ecef)

    # step 4: Calculate unit vector to satellites at previous time
    satellite_positions_previous = np.zeros((len(common_satellite_ids), 3))
    updated_pseudoranges_previous = np.zeros(len(common_satellite_ids))
    for j in range(len(common_satellite_ids)):
        svid = common_satellite_ids[j]
        old_pseudorange = measurements_previous[measurements_previous[:, 2] == svid][0, 3]
        sv_position, updated_pseudorange = calculate_sv_pos(previous_tor, old_pseudorange, svid, ephem_from_file)
        satellite_positions_previous[j] = sv_position
        updated_pseudoranges_previous[j] = updated_pseudorange
    user_position_previous_ecef, user_clock_error = calculate_least_squares_solution(updated_pseudoranges_previous, satellite_positions_previous)
    e_previous = np.zeros((len(common_satellite_ids), 3))
    for j in range(len(common_satellite_ids)):
        e_previous[j, :] = (satellite_positions_previous[j] - user_position_previous_ecef) / np.linalg.norm(satellite_positions_previous[j] - user_position_previous_ecef)

    # step 5: compute Doppler correction term
    svdoppler = np.zeros(len(common_satellite_ids))
    for j in range(len(svdoppler)):
        for dim in range(3):
            svdoppler += e_current[j, dim] * satellite_positions_current[j, dim] # TODO: should we really use these satellite positions?
            svdoppler -= e_previous[j, dim] * satellite_positions_previous[j, dim]

    # step 6: compute the Geometry correction term
    dgeom = np.zeros(len(common_satellite_ids))
    for j in range(len(dgeom)):
        for dim in range(3):
            svdoppler += e_previous[j, dim] * user_position_previous_ecef[dim]
            svdoppler -= e_current[j, dim] * user_position_previous_ecef[dim]

    # step 7: make Doppler- and Geometry adjustments
    z = np.zeros(len(common_satellite_ids))
    for j in range(len(z)):
        z[j] = delta_Phi[j] - svdoppler[j] - dgeom[j]

    # step 8: Set up H
    H = np.zeros((len(common_satellite_ids), 4))
    for j in range(len(H)):
        H[j, 0:3] = -e_current[j]
        H[j, 3] = 1

    # step 9: perform least squares
    delta_x = np.linalg.inv(H.T @ H) @ H.T @ z
    delta_x = delta_x / (current_tor - previous_tor)
    delta_x = delta_x[0:3]
    print(delta_x, np.linalg.norm(delta_x))
    velocities[time_counter, :] = delta_x
    times[time_counter] = current_tor
    time_counter += 1

velocities_lateral = np.zeros(len(velocities))
velocities_absolute = np.zeros(len(velocities))
for i in range(len(velocities)):
    velocities_lateral[i] = math.sqrt(velocities[i, 0]**2 + velocities[i, 1]**2)
    velocities_absolute[i] = math.sqrt(velocities[i, 0]**2 + velocities[i, 1]**2 + velocities[i, 2]**2)

plt.scatter(times, velocities_absolute)
plt.xlabel("time")
plt.ylabel("velocity (m/s)")
plt.title("Estimated velocity")
plt.show()


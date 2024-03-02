from matplotlib import pyplot as plt
from llh2ecef import *
from subroutines_for_position_estimation import *

raw = np.loadtxt("group04_raw.dat")
raw = raw[raw[:, 1] == 0] # only GPS data
raw = raw[raw[:, 5] == 0] # only L1 band
raw = raw[raw[:, 8] >= 3] # only data points with tracking status >= 3
ephem_from_file = np.loadtxt("group04_ephem.dat")

times = np.zeros(len(np.unique(raw[:, 0])))
user_positions_ecef = np.zeros((len(np.unique(raw[:, 0])), 3))
user_clock_errors = np.zeros(len(np.unique(raw[:, 0])))

for time_index in range(len(np.unique(raw[:, 0]))):
    tor = np.unique(raw[:, 0])[time_index]
    print("Time:", tor)
    data_points = raw[raw[:, 0] == tor]
    if len(data_points) < 4: continue # need at least four measurements to determine user position

    # now that we see at least four satellites, we can compute an actual user position
    sv_positions = np.zeros((len(data_points), 3))
    pseudoranges = np.zeros(len(data_points))
    for row in range(len(data_points)):
        sv_position, updated_pseudorange, svclock = calculate_sv_pos(tor, data_points[row, 3], data_points[row, 2], ephem_from_file)
        sv_positions[row] = sv_position
        pseudoranges[row] = updated_pseudorange
    user_position_ecef, user_clock_error = calculate_least_squares_solution(pseudoranges, sv_positions)

    times[time_index] = tor
    user_positions_ecef[time_index] = ecef2llh(user_position_ecef)
    user_clock_errors[time_index] = user_clock_error


plt.scatter(user_positions_ecef[:, 1], user_positions_ecef[:, 0], c=user_positions_ecef[:, 2])
plt.colorbar()
plt.xlabel("Longitude (degrees)")
plt.ylabel("Latitude (degrees)")
plt.title("Estimated positions in lat/lon/height (deg/deg/m)")

# plt.plot(times, user_positions_ecef[:,2])
# plt.ylabel("Height in m")
# plt.xlabel("Time of Reception")
# plt.title("Estimated height over time")

# plt.scatter(times, user_clock_errors)
# plt.ylabel("User clock offset [m]")
# plt.xlabel("Time of Reception")
# plt.title("Estimated error over time")
plt.show()


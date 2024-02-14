import numpy as np

raw = np.loadtxt("group01_raw.dat")
raw = raw[raw[:,1] == 0] # only GPS data
raw = raw[raw[:,5] == 0] # only L1 band
raw = raw[raw[:,8] >= 3] # only data points with tracking status >= 3

ephem = np.loadtxt("group01_ephem.dat")

print(raw)
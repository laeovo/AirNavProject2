import ephcal
import numpy as np
import math

ephem_data = np.loadtxt("ephem.dat")

for ephem_row in ephem_data:
    position, clock, error = ephcal()
    azimuth = math.atan2(position[0,0], position[1,0])
    elevation = math.atan2(position[2,0],math.sqrt(position[0,0]*position[0,0]+position[1,0]*position[1,0]))


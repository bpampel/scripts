#!/usr/bin/env python2
import glob
import numpy as np

# define filename and threshold for region of interest
kT = 2.49339  # 300K
F_max = 8*kT
outname = 'error_metric.txt'
error = []

# Find all fes files in folder (and sort them by time)
files = glob.glob('fes.b1.iter*')
time = [int(filter(str.isdigit, i)[1:]) for i in files]
files = [i for _, i in sorted(zip(time, files))]
time = sorted(time)



# Calculate reference by using the final FES
ref = np.genfromtxt(files[-1])
region = [ref[:, 1] < F_max]
F_values = np.extract(region, ref[:, 1])
ref[:, 1] -= np.average(F_values)

# perform analysis for all files (last one is to check if it works)
for i in range(len(files)):
    data = np.genfromtxt(files[i])
    F_values = np.extract(region, data[:, 1])
    data[:, 1] -= np.average(F_values)
    # calculate error metric using the reference FES
    variance = np.extract(region, (data[:, 1] - ref[:, -1])**2)
    error.append([time[i], np.sqrt(np.average(variance))])

np.savetxt(outname, np.asmatrix(error), fmt='%1.16f', delimiter=' ', newline='\n')

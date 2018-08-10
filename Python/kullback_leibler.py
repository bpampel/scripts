#!/usr/bin/env python2
import numpy as np
import glob
from scipy import stats

# define filename and threshold for region of interest
kT = 2.49339  # 300K
F_max = 8*kT
outname = 'kl_div.txt'
entropy = []

def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm

# Find all fes files in folder (and sort them by time)
files = glob.glob('fes.b1.iter*')
time = [int(filter(str.isdigit, i)[1:]) for i in files]
files = [i for _,i in sorted(zip(time, files))]
time = sorted(time)


# Retrieve reference distribution from the final FES
ref = np.genfromtxt(files[-1])
region = [ref[:,1] < F_max]
ref = normalize(np.extract(region, ref[:,1]))

# perform analysis for all files (last one is to check if it works)
for i in range(len(files)):
    data = np.genfromtxt(files[i])
    F_values = normalize(np.extract(region, data[:,1]))
    kldiv = stats.entropy(F_values, ref)
    # calculate error metric using the reference FES
    entropy.append([time[i], kldiv])

np.savetxt(outname, np.asmatrix(entropy), fmt='%1.16f', delimiter=' ', newline='\n')

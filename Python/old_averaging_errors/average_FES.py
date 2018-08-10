#!/usr/bin/env python3
"""
Averages pointwise over a set of FES and calculates the standart deviation.
Assumes the values of the sets are Gaussian distributed
"""
import glob
import os
import numpy as np

# get list of folders and fesfiles
folders = glob.glob('[0-9]*' + os.path.sep)
files = [os.path.basename(f) for f in glob.glob(folders[0]+'/fes.b1.iter*')]

if os.path.exists('avg'):
    os.rename('avg', 'bck.avg')
os.mkdir('avg')

reference = 'fes.ref.dat'

for filename in files:
    # set up array of right size from first folder
    tempdata = np.transpose(np.genfromtxt(folders[0]+filename))
        data = np.ndarray(shape=(len(folders), len(tempdata[0])), dtype=float)
    colvar = tempdata[0]
    data[0] = tempdata[1]

    # loop over all other folders
    for j in range(1, len(folders)):
        data[j] = np.transpose(np.genfromtxt(folders[j]+filename))[1]

    # all data is read, calculate averages and std dev
    avgdata = np.average(data, axis=0)
    stddev = np.std(data, axis=0, ddof=1, dtype=np.float64)

    # copy header from one infile
    comments = ''
    for line in open(folders[0]+filename):
        if line.startswith('#'):
            comments += line
        else:
            comments = comments[:-2] # remove last newline
            break

    # write data to file
    outdata = np.transpose(np.vstack((colvar, avgdata, stddev)))
    np.savetxt('avg/'+filename, np.asmatrix(outdata), header=comments,
               fmt='%1.16f', delimiter=' ', newline='\n')

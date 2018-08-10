#!/usr/bin/env python3
"""
Calculates mean and standart deviation for a set of retrieved FES
The averaging is performed at every point by bootstrapping
"""
import glob
import os
import time
from multiprocessing import Pool
import numpy as np


def calc_average(i):
    # set up arrays of right size from first folder
    tempdata = np.transpose(np.genfromtxt(folders[0] + files[i]))
    data = np.ndarray(shape=(len(folders), len(tempdata[0])), dtype=float)
    means = np.ndarray(shape=(len(tempdata[0])), dtype=float)
    stddev = np.ndarray(shape=(len(tempdata[0])), dtype=float)

    colvar = tempdata[0]
    data[0] = tempdata[1]

    # loop over all other folders
    for j in range(1, len(folders)):
        data[j] = np.transpose(np.genfromtxt(folders[j]+files[i]))[1]

    # all data is read, perform the calculations
    # random sampling for each CV value, calculate the means
    for k in range(0, len(colvar)):
        samples = np.array( [np.random.choice(data[:,k], len(data[:,k]))
                             for _ in range(5000)])
        means[k] = np.average(samples)
        stddev[k] = np.std(samples, ddof=1, dtype=np.float64)

    outdata = np.transpose(np.vstack((colvar, means, stddev)))
    np.savetxt('bootstr_avg/'+files[i], np.asmatrix(outdata),
               fmt='%1.16f', delimiter=' ', newline='\n')
    print('Averaging of FES file {} finished'.format(files[i]))
    return

if __name__ == '__main__':
    starttime = time.time()

    # get list of folders and fesfiles
    folders = glob.glob('[0-9]*' + os.path.sep)
    files = [os.path.basename(f) for f in glob.glob(folders[0]+'/fes.b1.iter*')]
    os.mkdir('bootstr_avg')
    try:
        os.system('taskset -p 0xff %d' % os.getpid())   # needed because of numpy
        pool = Pool(6)
        pool.map(calc_average, range(len(files)))
    finally:
        pool.close()
        pool.join()
        print('The job took {} seconds to complete'.format(time.time()-starttime))
else:
    exit(-1)

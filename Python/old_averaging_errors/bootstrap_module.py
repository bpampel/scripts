#!/usr/bin/env python2
"""
Calculates mean and standart deviation for a set of retrieved FES
The averaging is performed at every point by bootstrapping
"""
import glob
import os
import time
from multiprocessing import Pool
import numpy as np
import bootstrapped.bootstrap as bs
import bootstrapped.stats_functions as bs_stats


def calc_average(filename):
    # set up arrays of right size from first folder
    tempdata = np.transpose(np.genfromtxt(folders[0] + filename))
    data = np.ndarray(shape=(len(folders), len(tempdata[0])), dtype=float)
    means = np.ndarray(shape=(len(tempdata[0])), dtype=float)
    stddev = np.ndarray(shape=(len(tempdata[0])), dtype=float)

    colvar = tempdata[0]
    data[0] = tempdata[1]

    # loop over all other folders
    for i in range(1, len(folders)):
        data[i] = np.transpose(np.genfromtxt(folders[i]+filename))[1]

    # all data is read, perform the bootstrap
    for j in range(len(colvar)):
        means[j] = bs.bootstrap(data[:, j], stat_func=bs_stats.mean, num_iterations=5000).value
        stddev[j] = bs.bootstrap(data[:, j], stat_func=bs_stats.std, num_iterations=5000).value

    # copy header from one infile
    comments = ''
    for line in open(folders[0]+filename):
        if line.startswith('#'):
            comments += line
        else:
            comments = comments[:-2] # remove last newline
            break

    # write data to file
    outdata = np.transpose(np.vstack((colvar, means, stddev)))
    np.savetxt('bootstr_avg/'+filename, np.asmatrix(outdata), header=comments,
               fmt='%1.16f', delimiter=' ', newline='\n')
    print('Averaging of FES file {} finished'.format(filename))
    return

if __name__ == '__main__':
    starttime = time.time()

    # get list of folders and fesfiles
    folders = glob.glob('[0-9]*' + os.path.sep)
    files = [os.path.basename(f) for f in glob.glob(folders[0]+'/fes.b1.iter*')]

    if os.path.exists('bootstr_avg'):
        os.rename('bootstr_avg', 'bck.bootstr_avg')
    os.mkdir('bootstr_avg')

    try:
        os.system('taskset -p 0xff %d' % os.getpid())   # needed because of numpy
        pool = Pool(6)
        pool.map(calc_average, files)
    finally:
        pool.close()
        pool.join()
        print('The job took {} seconds to complete'.format(time.time()-starttime))
else:
    exit(-1)

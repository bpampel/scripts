#!/usr/bin/env python3
"""
Calculates error of a FES calculation
Averages pointwise over a set of FES and calculates the standart deviation.
Compares the average values pointwise to a reference FES for the bias.
Combination of both gives the total error.
"""

import glob
import os
import sys
import numpy as np


def get_filenames():
    """Returns all folders and fes files sorted by time"""
    tmp_folders = glob.glob('[0-9]*' + os.path.sep)
    tmp_files = [os.path.basename(f) for f in glob.glob(tmp_folders[0]+'/fes.b1.iter*')]
    tmp_times = [extract_time(f) for f in tmp_files]
    tmp_files = [i for _, i in sorted(zip(tmp_times, tmp_files))]
    tmp_times = sorted(tmp_times)
    return (tmp_folders, tmp_files, tmp_times)


def extract_time(x):
    """Returns time associated with file as int"""
    return int(''.join(i for i in x if i.isdigit())[1:])


def backup_if_exists(name):
    """Cascade of backups with format 'bck.$num.name'"""
    if os.path.exists(name):
        backupnum = 0
        while os.path.exists('bck.'+str(backupnum)+'.'+name):
            backupnum += 1
        os.rename(name, 'bck.'+str(backupnum)+'.'+name)


def extract_header(x):
    """Returns header of a plumed file as list"""
    header = []
    for line in open(x):
        if line.startswith('#'):
            header.append(line)
        else:
            return header


if __name__ == '__main__':
    # define some constants
    kT = 2.49339  # 300K
    shiftthreshold = 4*kT
    avgerror = []
    avgstddev = []
    avgbias = []

    if len(sys.argv) == 2:
        os.chdir(sys.argv[1])
    else:
        os.chdir(input("Base directory of the datasets: "))

    # read reference
    referencefile = 'fes.ref.dat'
    while True:
        try:
            ref = np.transpose(np.genfromtxt(referencefile))[1]
            break
        except IOError:
            referencefile = input("Path to the reference FES: ")

    # determine regions of interest
    shiftregion = [ref < shiftthreshold]
    errorregion = [ref < 2*shiftthreshold]
    refshift = np.average(np.extract(shiftregion, ref))

    folders, files, times = get_filenames()
    backup_if_exists('avg')
    os.mkdir('avg')

    for filename in files:
        # set up array of right size from first folder
        tempdata = np.transpose(np.genfromtxt(folders[0]+filename))
        data = np.ndarray(shape=(len(folders), len(tempdata[0])), dtype=float)
        colvar = tempdata[0]
        data[0] = tempdata[1]

        # loop over all other folders
        for j in range(1, len(folders)):
            data[j] = np.transpose(np.genfromtxt(folders[j]+filename))[1]

        # all data is read, shift data according to ref
        avgdata = np.average(data, axis=0)
        stddev = np.std(data, axis=0, ddof=1, dtype=np.float64)
        avgstddev.append(np.average(np.extract(errorregion, stddev)))

        datashift = np.average(np.extract(shiftregion, avgdata))
        avgdata = avgdata - datashift + refshift

        # calculate bias
        bias = np.absolute(avgdata - ref)
        avgbias.append(np.average(np.extract(errorregion, bias)))

        # add to recieve total error
        error = np.sqrt(stddev**2 + bias**2)
        avgerror.append(np.average(np.extract(errorregion, error)))

        # copy header from one infile and add fields
        fileheader = extract_header(folders[0]+filename)
        fileheader[0] = fileheader[0][:-1] + ' stddev bias error\n'
        # fileheader.append('#! shift ' + str(refshift + datashift) + '\n')
        fileheader = ''.join(fileheader)[:-1]

        # write data of current time to file
        outdata = np.transpose(np.vstack((colvar, avgdata, stddev, bias, error)))
        np.savetxt('avg/'+filename, np.asmatrix(outdata), header=fileheader,
                   comments='', fmt='%1.16f', delimiter=' ', newline='\n')

    # write averaged error to file
    backup_if_exists('error.txt')
    errorheader = '#! FIELDS time total_error stddev bias'
    errordata = np.transpose(np.vstack((times, avgerror, avgstddev, avgbias)))
    np.savetxt('error.txt', np.asmatrix(errordata), header=errorheader,
               comments='', fmt='%1.16f', delimiter=' ', newline='\n')

# else:
    # exit(-1)

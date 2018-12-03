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
    folders = glob.glob('[0-9]*' + os.path.sep)
    files = [os.path.basename(f) for f in glob.glob(folders[0]+'/fes.b1.iter*')]
    times = [extract_digits(f, 1) for f in files]
    files = [i for _, i in sorted(zip(times, files))]
    times = sorted(times)
    return (folders, files, times)


def extract_digits(x, start=0, end=None):
    """Returns digits contained in string as int.
    Arguments can be used if more than one number is contained
    This will result in an error if x does not contain any digit"""
    return int(''.join(i for i in x if i.isdigit())[start:end])


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

def get_nbins(header):
    """Extract the number of bins per direction from file header"""
    # watch out! If the corresponding CV label contains a number this fails
    return [extract_digits(s) for s in header if "nbin" in s]

def write_sliced_to_file(data, nbins, filename, header):
    """Writes the 2d data to file including a newline after every row"""
    data = data.reshape(*nbins, len(data[0])) # split into rows
    with open(filename, 'w') as outfile:
        outfile.write(header)
        for row in data[:-1]:
            np.savetxt(outfile, row, comments='', fmt='%1.16f',
                       delimiter=' ', newline='\n')
            outfile.write('\n')
        np.savetxt(outfile, data[-1], comments='', fmt='%1.16f',
                   delimiter=' ', newline='\n')


if __name__ == '__main__':
    # define some constants
    # kT = 2.49339  # 300K
    kT = 1
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
            ref = np.transpose(np.genfromtxt(referencefile))[2]
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
        tempdata = np.genfromtxt(folders[0]+filename).T
        data = np.ndarray(shape=(len(folders), len(tempdata[0])), dtype=float)
        colvar = tempdata[0], tempdata[1]
        data[0] = tempdata[2]

        # loop over all other folders
        for j in range(1, len(folders)):
            data[j] = np.transpose(np.genfromtxt(folders[j]+filename))[2]

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

        # copy header from one infile and check if bins match
        fileheader = extract_header(folders[0]+filename)
        nbins = get_nbins(fileheader)

        # modify header
        fileheader[0] = fileheader[0][:-1] + ' stddev bias error\n'
        fileheader = ''.join(fileheader)

        # write data of current time to file
        outdata = np.vstack((colvar, avgdata, stddev, bias, error)).T
        write_sliced_to_file(outdata, nbins, 'avg/'+filename, fileheader)

    # write averaged error to file
    backup_if_exists('error.txt')
    errorheader = '#! FIELDS time total_error stddev bias'
    errordata = np.vstack((times, avgerror, avgstddev, avgbias)).T
    np.savetxt('error.txt', np.asmatrix(errordata), header=errorheader,
               comments='', fmt='%1.16f', delimiter=' ', newline='\n')

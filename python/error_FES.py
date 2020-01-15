#!/usr/bin/env python3
"""
Calculates error of a FES calculation
Averages pointwise over a set of FES and calculates the standart deviation.
Compares the average values pointwise to a reference FES for the bias.
Combination of both gives the total error.
"""

import argparse
import glob
import os
import numpy as np


def parse_args():
    """Get cli args"""
    parser = argparse.ArgumentParser()
    parser.add_argument('path',
                        help="Path to the folder to be evaluated")
    parser.add_argument('-r', '--ref',
                        help="Path to the reference FES file",
                        required=True)
    parser.add_argument('-kT', '--kT', type=float,
                        help="Value of kT for the FES file (in matching units)",
                        required=True)
    parser.add_argument("-st", "--shifting_threshold", type=float,
                        help="Threshold value of FES (in units of kT) for shifting area",
                        default="4.0")
    parser.add_argument("-et", "--error_threshold", type=float,
                        help="Threshold value of FES (in units of kT) for error area",
                        default="8.0")
    # parser.add_argument("-min", "--minimum", type=float,
                        # help="Minimum of the CV range to be taken into account")
    # parser.add_argument("-max", "--maximum", type=float,
                        # help="Maximum of the CV range to be taken into account")
    args = parser.parse_args()
    return args


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
    args = parse_args()
    # define some constants and values
    shift_threshold = args.kT * args.shifting_threshold
    shift_threshold = args.kT * args.shifting_threshold
    error_threshold = args.kT * args.error_threshold
    # custom_cv_range = [0.0, 0.7]
    avgerror = []
    avgstddev = []
    avgbias = []
    avgnewbias = []
    avgstddevrms = []
    avgbiasrms = []
    avgnewbiasrms = []
    avgfolder = 'avg'
    errorfile = 'error.txt'

    # read reference
    try:
        colvar, ref = np.transpose(np.genfromtxt(args.ref))
    except IOError:
        print("Reference file not found!")

    try:
        os.chdir(args.path)
    except IOError:
        print("Could not change to specified directory")


    # determine regions of interest and create arrays of booleans
    # colvar_region = np.array([colvar >= custom_cv_range[0]]) & np.array([colvar <= custom_cv_range[1]])
    colvar_region = True # manual overwrite to include all
    shift_region = np.array([ref < shift_threshold]) & colvar_region
    error_region = np.array([ref < error_threshold]) & colvar_region
    errornorm = 1.0 / len(error_region)
    refshift = np.average(np.extract(shift_region, ref))

    folders, files, times = get_filenames()
    backup_if_exists(avgfolder)
    os.mkdir(avgfolder)

    for filename in files:
        # set up array of right size for all datasets
        data = np.ndarray(shape=(len(folders), len(colvar)), dtype=float)

        # loop over all folders
        for j in range(0, len(folders)):
            data[j] = np.transpose(np.genfromtxt(folders[j]+filename))[1]

        # all data is read, calculate error measures
        avgdata = np.average(data, axis=0)
        stddev = np.std(data, axis=0, ddof=1, dtype=np.float64)
        avgstddev.append(np.average(np.extract(error_region, stddev)))

        avgshift = np.average(np.extract(shift_region, avgdata))
        avgdata = avgdata - avgshift + refshift

        # calculate bias
        bias = np.absolute(avgdata - ref)
        avgbias.append(np.average(np.extract(error_region, bias)))

        # add to recieve total error
        error = np.sqrt(stddev**2 + bias**2)
        avgerror.append(np.average(np.extract(error_region, error)))

        # copy header from one infile and add fields
        fileheader = extract_header(folders[0]+filename)
        fileheader[0] = fileheader[0][:-1] + ' stddev bias error\n'
        # fileheader.append('#! shift ' + str(refshift + datashift) + '\n')
        fileheader = ''.join(fileheader)[:-1]

        # write data of current time to file
        outdata = np.transpose(np.vstack((colvar, avgdata, stddev, bias, error)))
        np.savetxt(avgfolder+ os.path.sep + filename, np.asmatrix(outdata), header=fileheader,
                   comments='', fmt='%1.16f', delimiter=' ', newline='\n')

    # write averaged error to file
    backup_if_exists(errorfile)
    errorheader =     "#! FIELDS time total_error stddev bias"
    errorheader += ("\n#! SET kT " + str(args.kT))
    errorheader += ("\n#! SET error_threshold " + str(args.error_threshold))
    errordata = np.transpose(np.vstack((times, avgerror, avgstddev, avgbias)))
    np.savetxt(errorfile, np.asmatrix(errordata), header=errorheader,
               comments='', fmt='%1.16f', delimiter=' ', newline='\n')

# else:
    # exit(-1)

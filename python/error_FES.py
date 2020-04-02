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
from helpers import misc as hlpmisc
from helpers import plumed_header as plmdheader


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
    parser.add_argument("--avgdir",
                        help='Name of the directory for the average file. Default is "avg"',
                        default="avg")
    parser.add_argument("--errorfile",
                        help='Name of the file for the error values. Default is "error.txt"',
                        default="error.txt")
    parser.add_argument("-st", "--shifting_threshold", type=float,
                        help="Threshold value of FES (in units of kT) for shifting area. Defaults to 4",
                        default="4.0")
    parser.add_argument("-et", "--error_threshold", type=float,
                        help="Threshold value of FES (in units of kT) for error area. Defaults to 8",
                        default="8.0")
    parser.add_argument('--cv-range', nargs='+', type=float,
                        help='CV range to be taken into account. Requires 2 values separated by spaces')
    args = parser.parse_args()
    if args.cv_range and len(args.cv_range) != 2:
        raise ValueError("--cv-range requires 2 values separated by spaces")
    return args



if __name__ == '__main__':
    args = parse_args()
    # define some constants and empty arrays for storage
    shift_threshold = args.kT * args.shifting_threshold
    shift_threshold = args.kT * args.shifting_threshold
    error_threshold = args.kT * args.error_threshold
    cv_region = True # full range by default
    avgerror = []
    avgstddev = []
    avgbias = []
    avgnewbias = []
    avgstddevrms = []
    avgbiasrms = []
    avgnewbiasrms = []


    # read reference
    try:
        colvar, ref = np.transpose(np.genfromtxt(args.ref))
    except IOError:
        print("Reference file not found!")


    # get folders and files
    folders = hlpmisc.get_subfolders(args.path)
    if len(folders) == 0:
        raise ValueError("No subfolders found at specified path.")
    files, times = hlpmisc.get_fesfiles(folders[0]) # assumes all folders have the same files


    # determine regions of interest and create arrays of booleans
    if args.cv_range:
        if args.cv_range[0] < colvar[0] or args.cv_range[1] > colvar[-1]:
            raise ValueError("Specified CV range is not contained in reference range [{}, {}]"
                             .format(colvar[0], colvar[-1]))
        cv_region = np.array([colvar >= args.cv_range[0]]) & np.array([colvar <= args.cv_range[1]])
    shift_region = np.array([ref < shift_threshold]) & cv_region
    error_region = np.array([ref < error_threshold]) & cv_region
    errornorm = 1.0 / len(error_region)
    refshift = np.average(np.extract(shift_region, ref))

    hlpmisc.backup_if_exists(args.avgfolder)
    os.mkdir(args.avgfolder)

    for filename in files:
        # set up array of right size for all datasets
        data = np.empty([len(folders), len(colvar)], dtype=float)

        # loop over all folders
        for j, folder in enumerate(folders):
            data[j] = np.transpose(np.genfromtxt(folder+filename))[1]

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
        fileheader = plmdheader.PlumedHeader()
        fileheader.parse_file(folders[0]+filename)
        fileheader[0] += ' stddev bias error'

        # write data of current time to file
        outdata = np.transpose(np.vstack((colvar, avgdata, stddev, bias, error)))
        np.savetxt(args.avgfolder+ os.path.sep + filename, np.asmatrix(outdata), header=fileheader,
                   comments='', fmt='%1.16f', delimiter=' ', newline='\n')

    # write averaged error to file
    hlpmisc.backup_if_exists(args.errorfile)
    errorheader =     "#! FIELDS time total_error stddev bias"
    errorheader += ("\n#! SET kT " + str(args.kT))
    errorheader += ("\n#! SET error_threshold " + str(args.error_threshold))
    errordata = np.transpose(np.vstack((times, avgerror, avgstddev, avgbias)))
    np.savetxt(args.errorfile, np.asmatrix(errordata), header=str(errorheader),
               comments='', fmt='%1.16f', delimiter=' ', newline='\n')

# else:
    # exit(-1)

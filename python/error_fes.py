#!/usr/bin/env python3
"""
Calculates error of a FES calculation
Averages pointwise over a set of FES and calculates the standart deviation.
Compares the average values pointwise to a reference FES for the bias.
Combination of both gives the total error.
"""

import argparse
from functools import partial
import os
from multiprocessing import Pool
import numpy as np
from helpers import misc as hlpmisc
from helpers import number_format as nfmt
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
    parser.add_argument("-o", "--outfile",
                        help='Name of the file for the error values. Default is "error.txt"',
                        default="error.txt")
    parser.add_argument("-st", "--shift-threshold", type=float, dest='shift_threshold',
                        help="Threshold value of FES (in units of kT) for shifting area. Defaults to 4",
                        default="4.0")
    parser.add_argument("-et", "--error-threshold", type=float, dest='error_threshold',
                        help="Threshold value of FES (in units of kT) for error area. Defaults to 8",
                        default="8.0")
    parser.add_argument('--cv-range', nargs='+', type=float,
                        help='CV range to be taken into account. Requires 2 values separated by spaces. \
                              Will be ignored for more than 2 dimensions')
    parser.add_argument("-np", "--numprocs", type=int, default="1",
                        help="Number of parallel processes")
    parser.add_argument("--individual-shifts", action='store_true', dest="individual_shifts",
                        help="Shift the individual FES instead of the average")
    args = parser.parse_args()
    if args.cv_range and len(args.cv_range) != 2:
        raise ValueError("--cv-range requires 2 values separated by spaces")
    return args


def calculate_error(filenames, avgdir, colvar, shift_region, error_region, ref, refshift, ind_shift):
    """
    Function doing the actual work.

    Arguments
    ---------
    filenames     : list with paths of all filenames to analyze
    avgdir        : directory to save averaged fes
    colvar        : numpy array holding the colvar values
    shift_region  : numpy array with booleans for the regions to consider for aligning data
    error_region  : numpy array with booleans for the regions to consider for error calculation
    ref           : numpy array holding the reference values
    ref_shift_avg : average value of ref in shift region

    Returns
    -------
    [avgstddev, avgbias, avgerror] : list of three floats
    avgstddev     : average standard deviation in error region
    avgbias       : average bias in error region
    avgerror      : combination of both error measures
    """
    dim, num_datapoints = colvar.shape
    data = np.empty([len(folders), num_datapoints], dtype=float)

    # read in all FES
    if ind_shift:
        for j, filename in enumerate(filenames):
            fes = np.transpose(np.genfromtxt(filename))[dim] # throw away colvar
            # exclude areas with inf in fes
            valid_fes_region = fes < np.inf
            valid_shift_region = np.bitwise_and(shift_region, valid_fes_region)
            refshift = np.average(ref[valid_shift_region])
            data[j] = fes - np.average(fes[valid_shift_region]) + refshift
    else:
        for j, filename in enumerate(filenames):
            data[j] = np.transpose(np.genfromtxt(filename))[dim] # throw away colvar

    # all data is read, calculate error measures
    avgdata = np.average(data, axis=0)
    stddev = np.std(data, axis=0, ddof=1, dtype=np.float64)
    avgstddev = np.average(np.extract(error_region, stddev))

    if not ind_shift:
        valid_avg_region = avgdata < np.inf
        valid_shift_region = np.bitwise_and(shift_region, valid_avg_region)
        refshift = np.average(ref[valid_shift_region])
        avgdata = avgdata - np.average(avgdata[valid_shift_region]) + refshift

    # calculate bias
    bias = np.absolute(avgdata - ref)
    avgbias = np.average(np.extract(error_region, bias))

    # add to recieve total error
    error = np.sqrt(stddev**2 + bias**2)
    avgerror = np.average(np.extract(error_region, error))

    # copy header from one infile and add fields
    fileheader = plmdheader.PlumedHeader()
    fileheader.parse_file(filenames[0])
    fileheader.fields += ["bias", "stddev", "total_error"]
    fileheader.set_constant("nruns_avg", len(filenames))
    fileheader.set_constant("shift_threshold", shift_threshold)
    fileheader.set_constant("error_threshold", error_threshold)

    # parse number format from FES file
    fmt_str = nfmt.get_string_from_file(filenames[0], dim)
    fmt = nfmt.NumberFmt(fmt_str)

    # write data of current time to file
    outfile = os.path.join(avgdir, os.path.basename(filenames[0]))
    outdata = np.transpose(np.vstack((colvar, avgdata, bias, stddev, error)))
    if dim == 1:
        np.savetxt(outfile, np.asmatrix(outdata), header=str(fileheader),
                   comments='', fmt=fmt.get(), delimiter=' ', newline='\n')
    if dim == 2:
        # find out number of bins per direction from header
        nbins = [val for key, val in fileheader.constants.items() if 'nbins' in key]
        write_sliced_to_file(outdata, nbins, outfile, fileheader, fmt.get())

    return [avgbias, avgstddev, avgerror]


def write_sliced_to_file(data, nbins, filename, header, fmts):
    """
    Writes 2d data to file including a newline after every row

    Arguments
    ---------
    data          : numpy array contaning positions and data information
    nbins         : list with bin numbers per direction
    filename      : path to write to
    header        : plumed_header to write to file
    fmts          : single format or list of formats for the data columns

    Returns
    -------
    Nothing
    """
    data = data.reshape(*nbins, len(data[0])) # split into rows
    with open(filename, 'w') as outfile:
        outfile.write(str(header) + '\n')
        for row in data[:-1]:
            np.savetxt(outfile, row, comments='', fmt=fmts,
                       delimiter=' ', newline='\n')
            outfile.write('\n')
        np.savetxt(outfile, data[-1], comments='', fmt=fmts,
                   delimiter=' ', newline='\n')


if __name__ == '__main__':
    args = parse_args()

    # define some constants and empty arrays for storage
    shift_threshold = args.kT * args.shift_threshold
    error_threshold = args.kT * args.error_threshold
    cv_region = True # full range by default
    fmt_times = '%10d'
    fmt_error = '%14.9f'

    # read reference and determine dimensions
    try:
        ref = np.genfromtxt(args.ref).T
    except IOError:
        print("Reference file not found!")
    dim = ref.shape[0] - 1
    colvar = ref[0:dim]
    ref = ref[dim]

    # get folders and files
    folders = hlpmisc.get_subfolders(args.path)
    if len(folders) == 0:
        raise ValueError("No subfolders found at specified path.")
    files, times = hlpmisc.get_fesfiles(folders[0]) # assumes all folders have the same files

    # determine regions of interest and create arrays of booleans
    if args.cv_range and dim==1: # missing implementation for higher dimensions
        if args.cv_range[0] < colvar[0] or args.cv_range[1] > colvar[-1]:
            raise ValueError("Specified CV range is not contained in reference range [{}, {}]"
                             .format(colvar[0], colvar[-1]))
        cv_region = np.bitwise_and(colvar >= args.cv_range[0], colvar <= args.cv_range[1])
    shift_region = np.bitwise_and(ref < shift_threshold, cv_region)
    error_region = np.bitwise_and(ref < error_threshold, cv_region)
    refshift = np.average(np.extract(shift_region, ref))

    avgdir = os.path.join(args.path, args.avgdir)
    hlpmisc.backup_if_exists(avgdir)
    os.mkdir(avgdir)

    # everything set up, now do the analysis for each time seperately
    filenames = [[os.path.join(d, f) for d in folders] for f in files]
    pool = Pool(processes=args.numprocs)
    avgvalues = pool.map(partial(calculate_error, avgdir=avgdir, colvar=colvar,
                                 shift_region=shift_region, error_region=error_region,
                                 ref=ref, refshift=refshift, ind_shift=args.individual_shifts),
                         filenames)

    # write averaged values to file
    errorfile = os.path.join(args.path, args.outfile)
    hlpmisc.backup_if_exists(errorfile)
    errordata = np.column_stack((times, avgvalues))
    fields = ["time", "bias", "stddev", "total_error"]
    errorheader = plmdheader.PlumedHeader(fields)
    errorheader.set_constant("kT", args.kT)
    errorheader.set_constant("shift_threshold", args.shift_threshold)
    errorheader.set_constant("error_threshold", args.error_threshold)
    fmt = [fmt_times] + 3*[fmt_error]
    np.savetxt(errorfile, np.asmatrix(errordata), header=str(errorheader),
               comments='', fmt=fmt, delimiter=' ', newline='\n')

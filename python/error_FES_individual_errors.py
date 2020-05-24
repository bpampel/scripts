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
from scipy.spatial.distance import pdist
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
    parser.add_argument('-kt', '--kT', type=float,
                        help="Value of kT for the FES file (in matching units)",
                        required=True)
    parser.add_argument("-o", "--outfile",
                        help='Name of the file for the error values. Default is "error.txt"',
                        default="error.txt")
    parser.add_argument("-st", "--shift-threshold", type=float, dest='shift_threshold',
                        help="Threshold value of FES (in units of kT) for shifting area. Defaults to 4",
                        default="4.0")
    parser.add_argument("-et", "--error-threshold", type=float, dest='error_threshold',
                        help="Threshold value of FES (in units of kT) for error area. Defaults to 8",
                        default="8.0")
    parser.add_argument('--cv-range', nargs='+', type=float, dest='cv_range',
                        help="CV range to be taken into account. Requires 2 values separated by spaces. \
                              Will be ignored for more than 2 dimensions")
    parser.add_argument("-np", "--numprocs", type=int, default="1",
                        help="Number of parallel processes")
    parser.add_argument("-em", "--error-metric", default="L1", dest='error_metric',
                        help="Distance metric for error calculation. \
                              Accepts 'L1' (manhattan) and 'L2' (euclidean). Defaults to 'L1'")
    args = parser.parse_args()
    if args.cv_range and len(args.cv_range) != 2:
        raise ValueError("--cv-range requires 2 values separated by spaces")
    return args


def calculate_error(filepath, dim, shift_region, error_region, ref, refshift, metric='L1'):
    """ Calculate error of FES wrt reference

    This shifts the read in FES by the average in the given region and then calculates the
    pointwise distance with a given metric. It returns the average distance in the region of interest.

    Arguments
    ---------
    filepath      : paths to FES to analyze
    dim           : dimensions of FES
    shift_region  : numpy array with booleans for the regions to consider for aligning data
    error_region  : numpy array with booleans for the regions to consider for error calculation
    ref           : numpy array holding the reference values
    refshift      : average value of ref in shift region
    metric        : used metric for distance calculation. Accepts either 'L1' or 'L2'.

    Returns
    -------
    error         : float with the average error value in error_region
    """
    fes = np.transpose(np.genfromtxt(filepath))[dim]  # throw away colvar
    fes = fes - np.average(np.extract(shift_region, fes)) + refshift
    if metric == 'L1':
        error = np.abs(fes - ref)
    elif metric == 'L2':
        error = np.sqrt((fes - ref)**2)
    else:
        raise ValueError("Metric must be 'L1' or 'L2'")
    return np.average(np.extract(error_region, error))


def main():
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
    files, times = hlpmisc.get_fesfiles(folders[0])  # assumes all folders have the same files

    # determine regions of interest and create arrays of booleans
    if args.cv_range and dim==1: # missing implementation for higher dimensions
        if args.cv_range[0] < colvar[0] or args.cv_range[1] > colvar[-1]:
            raise ValueError("Specified CV range is not contained in reference range [{}, {}]"
                             .format(colvar[0], colvar[-1]))
        cv_region = np.array([colvar >= args.cv_range[0]]) & np.array([colvar <= args.cv_range[1]])
    shift_region = np.array([ref < shift_threshold]) & cv_region
    error_region = np.array([ref < error_threshold]) & cv_region
    refshift = np.average(np.extract(shift_region, ref))

    # everything set up, now calculate errors for all files
    filepaths = [os.path.join(d, f) for d in folders for f in files]
    pool = Pool(processes=args.numprocs)
    errors = pool.map(partial(calculate_error, dim=dim, shift_region=shift_region,
                              error_region=error_region, ref=ref, refshift=refshift,
                              metric=args.error_metric), filepaths)
    errors = np.array(errors).reshape(len(folders),len(files))  # put in matrix form

    # write error for each folder to file
    fileheader = plmdheader.PlumedHeader()
    fileheader.add_line("FIELDS time error")
    fileheader.add_line("SET kT {}".format(args.kT))
    fileheader.add_line("SET shift_threshold {}".format(args.shift_threshold))
    fileheader.add_line("SET error_threshold {}".format(args.error_threshold))
    fmt = [fmt_times] + [fmt_error]
    for i, folder in enumerate(folders):
        errorfile = os.path.join(folder, args.outfile)
        hlpmisc.backup_if_exists(errorfile)
        np.savetxt(errorfile, np.vstack((times, errors[i])).T, header=str(fileheader),
                   comments='', fmt=fmt, delimiter=' ', newline='\n')

    # calculate average and stddev
    avg_error = np.average(errors, axis=0)
    stddev = np.std(errors, axis=0, ddof=1)
    # write to file
    avgfile = os.path.join(args.path, args.outfile)  # in base dir
    hlpmisc.backup_if_exists(avgfile)
    fileheader[0] = "FIELDS time avg_error stddev"
    fileheader.add_line('SET nruns_avg {}'.format(len(filepaths)))
    fmt.append(fmt_error)
    np.savetxt(avgfile, np.vstack((times, avg_error, stddev)).T, header=str(fileheader),
               comments='', fmt=fmt, delimiter=' ', newline='\n')


if __name__ == '__main__':
    main()
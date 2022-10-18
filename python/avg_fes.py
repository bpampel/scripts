#!/usr/bin/env python3
"""
Calculates error of a FES calculation
Averages pointwise over a set of FES and calculates the standart deviation.
Compares the average values pointwise to a reference FES for the bias.
Combination of both gives the total error.
"""

import argparse
import numpy as np
import sys
from helpers import misc as hlpmisc
from helpers import number_format as nfmt
from helpers import plumed_header as plmdheader


def parse_args():
    """Get cli args"""
    parser = argparse.ArgumentParser()
    parser.add_argument("files",
                        help="Path to the files to be evaluated", nargs="+")
    parser.add_argument("-d", "--dim", type=int, default="1",
                        help="Dimensions of FES, defaults to 1")
    parser.add_argument("-kT", "--kT", type=float,
                        help="Value of kT for the FES file (in matching units)",
                        required=True)
    parser.add_argument("-o", "--outfile",
                        help='Name of the file for the averaged fes. Default is "fes_avg"',
                        default="fes_avg")
    parser.add_argument("-st", "--shifting_threshold", type=float,
                        help="Threshold value of FES (in units of kT) for shifting area. Defaults to 4",
                        default="4.0")
    parser.add_argument("-np", "--numprocs", type=int, default="1",
                        help="Number of parallel processes")
    return parser.parse_args()


def avg_fes(filenames, outfile, colvar, shift_region, refshift):
    """Calculate the average FES and save it to file

    Arguments
    ---------
    filenames     : list with paths of all filenames to analyze
    outfile       : path for averaged fes
    colvar        : numpy array holding the colvar values
    shift_region  : numpy array with booleans for the regions to consider for aligning data
    refshift      : average value of ref in shift region

    Returns
    -------
    Nothing
    """
    dim, num_datapoints = colvar.shape
    data = np.empty([len(filenames), num_datapoints], dtype=float)

    # loop over all folders
    for j, filename in enumerate(filenames):
        tmp = np.transpose(np.genfromtxt(filename))[dim] # throw away colvar
        data[j] = tmp - np.average(np.extract(shift_region, tmp)) + refshift

    avgdata = np.average(data, axis=0)

    # copy header from one infile and add fields
    fileheader = plmdheader.create_from_file(filenames[0])
    fileheader.set_constant("nruns_avg", len(filenames))

    # parse number format from FES file
    fmt_str = nfmt.get_string_from_file(filenames[0], dim)
    fmt = nfmt.NumberFmt(fmt_str)

    # write data of current time to file
    outdata = np.transpose(np.vstack((colvar, avgdata)))
    if dim == 1:
        np.savetxt(outfile, np.asmatrix(outdata), header=str(fileheader),
                   comments='', fmt=fmt.get(), delimiter=' ', newline='\n')
    if dim == 2:
        # find out number of bins per direction from header constants
        nbins = [int(val) for key, val in fileheader.constants.items() if 'nbins' in key]
        hlpmisc.write_2d_sliced_to_file(outfile, outdata, nbins, fmt.get(), fileheader)


def main():
    args = parse_args()

    # define some constants
    shift_threshold = args.kT * args.shifting_threshold

    n_files=len(args.files)
    if n_files <= 1:
        raise ValueError("You must specify at least 2 FES files to average")
        sys.exit(-1)

    # use first as reference:
    # determine dimensions and shifting
    try:
        fes = np.genfromtxt(args.files[0]).T
    except IOError:
        print(f"File {args.files[0]} could not be read!")
        sys.exit(-1)
    colvar = fes[0:args.dim]
    ref = fes[args.dim]
    shift_region = np.array([ref < shift_threshold])
    refshift = np.average(np.extract(shift_region, ref))

    # now average, it is assumed that all files have exactly the same structure as the first
    avg_fes(args.files, args.outfile, colvar, shift_region, refshift)


if __name__ == '__main__':
    main()

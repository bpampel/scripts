#!/usr/bin/env python3
"""
Project a multidimensional free energy surface onto fewer dimensions.
Currently only 2d -> 1d
"""

import argparse
from os import path
import numpy as np
from helpers import number_format as nfmt


def parse_args():
    """Get cli args"""
    parser = argparse.ArgumentParser()
    parser.add_argument('filename',
                        help="Path to the FES file to be projected")
    parser.add_argument("-d", "--dim", type=int,
                        help="Dimension on which should be projected\
                              (given as colum number of the FES file)",
                        required=True)
    parser.add_argument("-kT", "--kT", type=float,
                        help="Energy (in units of kT) of the FES file",
                        required=True)
    parser.add_argument("-o", "--outfile",
                        help="Name of the output file")

    args = parser.parse_args()

    if args.outfile is None:
        directory = ''
        splitted_filename = args.filename.split('/')
        if len(splitted_filename) > 1:
            directory = '/'.join(splitted_filename[:-1]) + "/"
        args.outfile = directory + "proj_" + str(args.dim) + "_" +  splitted_filename[-1]

    return args






# def extract_header(x):
    # """Returns header of a plumed file as list"""
    # header = []
    # for line in open(x):
        # if line.startswith('#'):
            # header.append(line)
        # else:
            # return header

def get_number_string(filename):
    """Get number string of file from the first number line"""
    for line in open(filename):
        if line.startswith('#'):
            continue
        else:
            return line.split()[2]


if __name__ == '__main__':
    # define some constants and values
    args = parse_args()

    fes = np.genfromtxt(args.filename)
    fmt = nfmt.NumberFmt(get_number_string(args.filename))

    num_dim1 = np.where(fes[:, 0] == fes[0, 0])[0][1] # via periodicity
    num_dim2 = int(fes.shape[0] / num_dim1)

    # get CV values of correct dimension
    if args.dim == 1:
        cv_values = fes[:num_dim1, 0]
    elif args.dim == 2:
        cv_values = fes[::num_dim1, 1]
    else:
        raise ValueError("Dimension must be either 1 or 2")

    # put values in matrix, dim2 is in the rows (because fes from plumed is column major)
    fesmatrix = fes[:, 2].reshape(num_dim2, num_dim1)
    probabilities = np.exp(- fesmatrix / float(args.kT))

    cv_delta = cv_values[1] - cv_values[0]
    projected_probabilities = np.trapz(probabilities, axis=args.dim-1, dx=cv_delta)
    projected_fes = - args.kT * np.log(projected_probabilities)
    projected_fes -= np.min(projected_fes)

    np.savetxt(args.outfile, np.asmatrix(np.vstack((cv_values, projected_fes)).T),
               fmt=fmt.get(), delimiter=' ', newline='\n')

# else:
    # exit(-1)

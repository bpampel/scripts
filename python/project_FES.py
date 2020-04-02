#!/usr/bin/env python3
"""
Project a multidimensional free energy surface onto fewer dimensions.
Currently only 2d -> 1d
"""

import argparse
import numpy as np
from helpers import plumed_header as plmdheader
from helpers import number_format as nfmt
from helpers import misc as hlpmisc


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
        args.outfile = hlpmisc.prefix_filename(args.filename, "proj_" + str(args.dim) + "_")

    return args


def manipulate_header(header, dim):
    """
    Change the original header of the input file
    Remove everything from the projected out dimension
    """
    fieldslinenum, fields = header.search_lines("FIELDS")[0]
    proj_variable = fields.split()[dim+1]
    proj_value = fields.split()[4]
    header[fieldslinenum] = "#! FIELDS {} projection.{}".format(proj_variable, proj_value)
    removed_value = fields.split()[4-dim]
    header.del_lines([i for i, _ in header.search_lines(removed_value)])


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
    header = plmdheader.PlumedHeader()
    header.parse_file(args.filename)
    manipulate_header(header, args.dim)


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
               header=str(header), fmt=fmt.get(), comments='', delimiter=' ',
               newline='\n')

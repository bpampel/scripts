#!/usr/bin/env python3
"""
Project a multidimensional free energy surface onto fewer dimensions.
Currently only 2d -> 1d
"""

import argparse
from os import path
import numpy as np


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


def number_format_parser(number):
    """Get the format of a number string"""
    # all variables here, could be turned into module
    fmt_char = 'f'
    is_negative = False
    int_length = 0
    dec_length = 0

    try:
        float(number)
    except ValueError:
        print('String is not a number')

    if number[0] == '-':
        is_negative = True

    dot_position = number.find('.')
    if dot_position == -1:
        fmt_str = "%{}.0f".format(len(number) - is_negative)
        return fmt_str

    int_length = dot_position - is_negative

    # check if exponential if dot was in 2nd (with minus 3rd) character
    if dot_position - is_negative == 1:
        for e in ['e', 'E']:
            e_position = number.find(e)
            if e_position != -1:
                fmt_char = e
                dec_length = e_position - dot_position - 1
                break

    if fmt_char == 'f':
        dec_length = len(number) - dot_position - 1

    fmt_str = "%{}.{}{}".format(int_length, dec_length, fmt_char)
    return fmt_str




# def extract_header(x):
    # """Returns header of a plumed file as list"""
    # header = []
    # for line in open(x):
        # if line.startswith('#'):
            # header.append(line)
        # else:
            # return header

def get_number_format(filename):
    """Get number format of file from the first number"""
    for line in open(filename):
        if line.startswith('#'):
            continue
        else:
            return number_format_parser(line.split()[2])


if __name__ == '__main__':
    # define some constants and values
    args = parse_args()

    fes = np.genfromtxt(args.filename)
    number_format = get_number_format(args.filename)

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
               fmt=number_format, delimiter=' ', newline='\n')

# else:
    # exit(-1)

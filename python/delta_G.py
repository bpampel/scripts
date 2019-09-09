#!/usr/bin/env python3
"""
Calculate free energy difference of two states by integrating over the
probabilities
"""

import argparse
import numpy as np
from helpers import misc as hlpmisc


def parse_args():
    """Get cli args"""
    parser = argparse.ArgumentParser()
    parser.add_argument('path',
                        help="Path to the FES file or folder to be evaluated")
    parser.add_argument('-kT', '--kT', type=float,
                        help="Energy (in units of kT) of the FES file",
                        required=True)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-f', '--file', action='store_false',
                       help="Parse single fes file (default)")
    group.add_argument('-d', '--dir', action='store_true',
                       help="Parse one or multiple directories. \
                             Looks for all numbered [0-9]* subfolders \
                             or works on the directory itself.")
    # parser.add_argument("-A", "--stateA", '-nargs', nargs='+', type=float,
                        # help="Approximate location of basin A (takes 2 values)")
    # parser.add_argument("-B", "--stateB", '-nargs', nargs='+', type=float,
                        # help="Approximate location of basin B (takes 2 values)")
    parser.add_argument("-t", "--threshold", type=float,
                        help="Probability threshold of basins",
                        default="0.0")
    parser.add_argument("-o", "--outfile",
                        help="Name of the output file")

    args = parser.parse_args()

    if args.outfile is None:
        directory = ''
        splitted_filename = args.filename.split('/')
        if len(splitted_filename) > 1:
            directory = '/'.join(splitted_filename[:-1]) + "/"
        args.outfile = directory + "delta_G" + "_" +  splitted_filename[-1]

    return args


def convert_to_matrix(fes):
    num_dim1 = np.where(fes[:, 0] == fes[0, 0])[0][1] # via periodicity
    num_dim2 = int(fes.shape[0] / num_dim1)

    # dim2 is in the rows because fes from plumed is column major
    return fes[:, 2].reshape(num_dim2, num_dim1)


def calculate_delta_G(filename, kT, masks):
    """
    Calculates the free energy difference between two states

    Arguments
    ---------
    filename : path to the fes file
    kT       : energy in units of kT
    masks    : a list with boolean numpy arrays resembling the two states

    Returns
    -------
    delta_G  : a double containing the free energy difference
    """

    fes = np.genfromtxt(filename)
    fesmatrix = convert_to_matrix(fes)

    probabilities = np.exp(- fesmatrix / float(kT))

    # masks for upper and lower regions of CV space
    probs = []

    masks = masks & (probabilities > args.threshold) # exclude low probability areas

    for i in range(2):
        probs.append(np.sum(probabilities[masks[i]]))

    delta_G = - kT * np.log(probs[0]/probs[1])

    return(delta_G)


if __name__ == '__main__':
    # read in cli arguments
    args = parse_args()

    # could also read in masks from elsewhere, this is for the Wulfe-Quapp potential
    masks = []
    masks.append(np.vstack((np.full((150, 301), True), np.full((151, 301), False))))
    masks.append(np.vstack((np.full((151, 301), False), np.full((150, 301), True))))

    if args.file:
        print(calculate_delta_G(args.path, args.kT, masks))
    elif args.dir:
        folders = hlpmisc.get_subfolders(args.path)
        if not folders: # no subdirectories found - use current one
            folders = [args.path]
        files, times = hlpmisc.get_fesfiles(folders[0])
        for folder in folders:
            for filename, i in enumerate(files):
                print(calculate_delta_G(filename, args.kT, masks))



# else:
    # exit(-1)

#!/usr/bin/env python3
"""
Calculate free energy difference of two states by integrating over the
probabilities
"""

import argparse
import numpy as np
from helpers.misc import get_filenames


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


if __name__ == '__main__':
    # define some constants and values
    args = parse_args()

    fes = np.genfromtxt(args.filename)
    fesmatrix = convert_to_matrix(fes)

    cv1 = fes[:fesmatrix.shape[0], 0]
    cv2 = fes[::fesmatrix.shape[0], 1]

    probabilities = np.exp(- fesmatrix / float(args.kT))

    # masks for upper and lower regions of CV space
    masks = []
    probs = []
    masks.append(np.vstack((np.full((150, 301), True), np.full((151, 301), False))))
    masks.append(np.vstack((np.full((151, 301), False), np.full((150, 301), True))))

    masks = masks & (probabilities > args.threshold) # exclude low probability areas

    for i in range(2):
        probs.append(np.sum(probabilities[masks[i]]))

    delta_G = - args.kT * np.log(probs[0]/probs[1])

    print(delta_G)




# else:
    # exit(-1)

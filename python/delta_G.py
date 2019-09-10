#!/usr/bin/env python3
"""
Calculate free energy difference of two states by integrating over the
probabilities
"""

import argparse
import numpy as np
import os
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
    group.add_argument('-f', '--file', action='store_const', dest='fd', const='f',
                       help="Parse single fes file (default)")
    group.add_argument('-d', '--dir', action='store_const', dest='fd', const='d',
                       help="Parse one or multiple directories. \
                             Looks for all numbered [0-9]* subfolders \
                             or works on the directory itself.")
    parser.add_argument("-avg", "--average",
                        help="Also calculate average over runs. \
                              Requires the --dir flag.")
    # parser.add_argument("-A", "--stateA", '-nargs', nargs='+', type=float,
                        # help="Approximate location of basin A (takes 2 values)")
    # parser.add_argument("-B", "--stateB", '-nargs', nargs='+', type=float,
                        # help="Approximate location of basin B (takes 2 values)")
    parser.add_argument("-t", "--threshold", type=float,
                        help="Probability threshold of basins",
                        default="0.0")
    parser.add_argument("-o", "--outfile",
                        help='Name of the output file(s). Default is "delta_G"')
    parser.set_defaults(fd='f')

    args = parser.parse_args()

    if args.average and not args.dir:
        raise ValueError("-avg without -d doesn't make sense.")

    if args.fd == 'f' and args.outfile:
        print("Single file given. Ignoring outfile argument")


    return args


def get_outfilenames(outfile, folders, inputpath):
    """decide on outfile names depending on cli arguments given"""
    outfilenames = []
    avgname = None
    if outfile:
        if len(folders) == 1:
            # next two lines would put the file by default in the input folder. Desired behaviour?
            # if not os.path.dirname(outfile):
                # outfile = os.path.dirname(inputpath) + outfile
            outfilenames.append(outfile)
        else:
            dirname = os.path.dirname(outfile)
            if not dirname: # if only filename without path
                outfilenames = [folder + outfile for folder in folders]
                avgname = os.path.join(*dirname.split(os.path.sep)[:-1], outfile) # same name in base directory
            else: # put all files in given folder with numbers to differentiate
                dirname = os.path.dirname(outfile)
                basename = os.path.basename(outfile)
                numbers = [folder.split(os.path.sep)[-2] for folder in folders]
                outfilenames = [os.path.join(dirname, basename + "." + i) for i in numbers]
                avgname = dirname + os.path.sep + basename + ".avg"
    else: # no filename given
        basename = "delta_G"
        outfilenames = [folder + basename for folder in folders]
        avgname = os.path.join(*dirname.split(os.path.sep)[:-1], "delta_G") # same name in base directory


    return outfilenames


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

    # masks = masks & (probabilities > args.threshold) # exclude low probability areas

    for i in range(2):
        probs.append(np.sum(probabilities[masks[i]]))

    delta_G = - kT * np.log(probs[0]/probs[1])

    return delta_G


def main():
    # read in cli arguments
    args = parse_args()


    # could also read in masks from elsewhere, this is for the Wulfe-Quapp potential
    masks = []
    masks.append(np.vstack((np.full((150, 301), True), np.full((151, 301), False))))
    masks.append(np.vstack((np.full((151, 301), False), np.full((150, 301), True))))

    if args.fd == 'f':
        print(calculate_delta_G(args.path, args.kT, masks))

    elif args.fd == 'd':
        folders = hlpmisc.get_subfolders(args.path)

        if not folders: # no subdirectories found - use only given one
            if args.path[-1] != os.path.sep:
                args.path += os.path.sep # add possibly missing "/"
            folders = [args.path]
            if args.average:
                raise ValueError("No subdirectories found. Averaging not possible.")

        # has to be done here as it's previously not clear if multiple folders are involved
        outfilenames, avgfilename = get_outfilenames(args.outfile, folders, args.path)

        files, times = hlpmisc.get_fesfiles(folders[0])

        delta_G = np.ndarray(shape=(len(folders),len(files))) # store all for stddev

        errorheader = '#! FIELDS time delta_G stddev'

        for i, folder in enumerate(folders):
            for j, filename in enumerate(files):
                delta_G[i,j] = calculate_delta_G(folder + filename, args.kT, masks)
            np.savetxt(outfilenames[i], delta_G[i])

        if args.average:
            avg_delta_G = np.average(delta_G, axis=0)
            stddev = np.std(delta_G, axis=0, ddof=1)
            np.savetxt(avgfilename, np.vstack(times, avg_delta_G, stddev))


if __name__ == '__main__':
    main()

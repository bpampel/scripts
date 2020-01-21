#!/usr/bin/env python3
"""
Calculate free energy difference of two states by integrating over the
probabilities
"""

import argparse
from functools import partial
import os
from multiprocessing import Pool
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
    group.add_argument('-f', '--file', action='store_const', dest='fd', const='f',
                       help="Parse single fes file (default)")
    group.add_argument('-d', '--dir', action='store_const', dest='fd', const='d',
                       help="Parse one or multiple directories. \
                             Looks for all numbered [0-9]* subfolders \
                             or works on the directory itself.")
    parser.add_argument("-avg", "--average", action='store_true',
                        help="Also calculate average over runs. \
                              Requires the --dir flag.")
    parser.add_argument("-np", "--numprocs", type=int, default="1",
                        help="Number of parallel processes")
    parser.add_argument("-o", "--outfile",
                        help='Name of the output file(s). Default is "delta_G"')
    parser.set_defaults(fd='f')

    args = parser.parse_args()

    # exclude some options that don't make sense together
    if args.fd == 'f':
        if args.average:
            raise ValueError("-avg without -d doesn't make sense.")
        if args.outfile:
            print("Single file given. Ignoring outfile argument.")
    return args


def get_outfilenames(outfile, folders):
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
                avgname = os.path.join(os.path.dirname(os.path.dirname(folders[0])), outfile) # same name in base directory
            else: # put all files in given folder with numbers to differentiate
                dirname = os.path.dirname(outfile)
                basename = os.path.basename(outfile)
                numbers = [folder.split(os.path.sep)[-2] for folder in folders]
                outfilenames = [os.path.join(dirname, basename + "." + i) for i in numbers]
                avgname = os.path.join(dirname, basename + ".avg")
    else: # no filename given
        basename = "delta_G"
        outfilenames = [folder + basename for folder in folders]
        avgname = os.path.join(os.path.dirname(os.path.dirname(folders[0])), "delta_G") # same name up one directory


    return (outfilenames, avgname)



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

    fes = np.genfromtxt(filename).T[1]

    probabilities = np.exp(- fes / float(kT))

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


    # could also read in masks from elsewhere, this is for the Wolfe-Quapp potential
    masks = []
    masks.append(np.concatenate((np.full(150, True), np.full(151, False))))
    masks.append(np.concatenate((np.full(151, False), np.full(150, True))))

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
        outfilenames, avgfilename = get_outfilenames(args.outfile, folders)

        files, times = hlpmisc.get_fesfiles(folders[0])

        pool = Pool(processes=args.numprocs)

        allfilenames = [folder + file for folder in folders for file in files]
        delta_G = pool.map(partial(calculate_delta_G, kT=args.kT, masks=masks), allfilenames)
        delta_G = np.array(delta_G).reshape(len(folders),len(files))

        for i, folder in enumerate(folders):
            np.savetxt(outfilenames[i], np.vstack((times, delta_G[i])).T)

        if args.average:
            avg_delta_G = np.average(delta_G, axis=0)
            stddev = np.std(delta_G, axis=0, ddof=1)
            np.savetxt(avgfilename, np.vstack((times, avg_delta_G, stddev)).T)

if __name__ == '__main__':
    main()

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
from helpers import plumed_header as plmdheader


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
    parser.add_argument("-m", "--masks", type=argparse.FileType('r'), nargs='+',
                        help="Files containing masks of first basin")
    parser.add_argument("-np", "--numprocs", type=int, default="1",
                        help="Number of parallel processes")
    parser.add_argument("-o", "--outfile",
                        help='Name of the output file(s). Default is "delta_F"')
    parser.set_defaults(fd='f')

    args = parser.parse_args()

    # exclude some options that don't make sense together
    if args.fd == 'f' and args.outfile:
        print("Single file given. Ignoring outfile argument.")
    return args


def get_outfilenames(outfile, folders):
    """
    decide on outfile names depending on cli arguments given

    Arguments
    ---------
    outfile  : user specified outfile name/path
    folders  : list with all folders to evaluate

    Returns
    -------
    (outfilenames, avgname) : tuple
    outfilenames            : list containing the outfile name for each folder
    avgname                 : string with the path for the avg file
    """
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
            if not dirname:  # if only filename without path
                outfilenames = [os.path.join(d, outfile) for d in folders]
                avgname = os.path.join(os.path.dirname(os.path.dirname(folders[0])), outfile) # same name in base directory
            else:  # put all files in given folder with numbers to differentiate
                basename = os.path.basename(outfile)
                numbers = [folder.split(os.path.sep)[-2] for folder in folders]
                outfilenames = [os.path.join(dirname, basename + "." + i) for i in numbers]
                avgname = os.path.join(dirname, basename + ".avg")
    else:  # no filename given
        basename = "delta_F"
        outfilenames = [folder + basename for folder in folders]
        avgname = os.path.join(os.path.dirname(os.path.dirname(folders[0])), "delta_F")  # same name up one directory

    return (outfilenames, avgname)


def calculate_state_probabilities(filename, kT, masks):
    """
    Calculates the free energy difference between two states

    Arguments
    ---------
    filename : path to the fes file
    kT       : energy in units of kT
    masks    : a list of numpy arrays resembling the states

    Returns
    -------
    probs    : a list of doubles containing the probabilities of the states
    """

    fes = np.genfromtxt(filename).T[-1]  # assumes that last column is free energy
    if len(fes) != len(masks[0]):
        raise ValueError('Masks and FES of file {} are not of the same length ({} and {})'
                         .format(filename, len(masks[0]), len(fes)))

    probabilities = np.exp(- fes / float(kT))
    norm = np.sum(probabilities)
    state_probs = [np.sum(probabilities[m])/norm for m in masks]
    return state_probs


def main():
    # read in cli arguments, define constants
    args = parse_args()
    fmt_times = '%10d'
    fmt_probs = '%14.9f'

    masks = []
    for m in args.masks:
        try:
            mask = np.genfromtxt(m).astype('bool')  # could also save in binary but as int/bool is more readable
        except OSError:
            print('Error: Specified masks file "{}" not found'.format(m))
            raise
        masks.append(mask)
    if not all(len(m) == len(masks[0]) for m in masks[1:]):
        raise ValueError('Not all masks are of the same length! Length are {}'
                         .format([len(m) for m in masks]))

    if args.fd == 'f':
        probs = calculate_state_probabilities(args.path, args.kT, masks)
        for p in probs:
            print(fmt_probs % p)

    elif args.fd == 'd':
        folders = hlpmisc.get_subfolders(args.path)

        if not folders:  # no subdirectories found - use only given one
            if args.path[-1] != os.path.sep:
                args.path += os.path.sep  # add possibly missing "/"
            folders = [args.path]

        # has to be done here as it's previously not clear if multiple folders are involved
        outfilenames, avgfilename = get_outfilenames(args.outfile, folders)

        files, times = hlpmisc.get_fesfiles(folders[0])

        pool = Pool(processes=args.numprocs)

        allfilenames = [os.path.join(d, f) for d in folders for f in files]
        probs = pool.map(partial(calculate_state_probabilities, kT=args.kT, masks=masks), allfilenames)
        probs = np.array(probs).reshape(len(folders), len(files), len(masks))

        header = plmdheader.PlumedHeader()
        fields = 'FIELDS time'
        for i,_ in enumerate(masks):
            fields += ' mask' + str(i)
        header.add_line(fields)
        header.add_line('SET kT {}'.format(args.kT))
        fmt = [fmt_times] + [fmt_probs] * len(masks)

        for i, f in enumerate(outfilenames):
            hlpmisc.backup_if_exists(f)
            np.savetxt(f, np.vstack((times, probs[i].T)).T, header=str(header), fmt=fmt,
                       comments='', delimiter=' ', newline='\n')

if __name__ == '__main__':
    main()

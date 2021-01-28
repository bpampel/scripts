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
    parser.add_argument('-avg', '--average', action='store_true',
                        help="Also calculate average over runs. \
                              Requires the --dir flag.")
    # parser.add_argument("-A", "--stateA", '-nargs', nargs='+', type=float,
                        # help="Approximate location of basin A (takes 2 values)")
    # parser.add_argument("-B", "--stateB", '-nargs', nargs='+', type=float,
                        # help="Approximate location of basin B (takes 2 values)")
    # parser.add_argument("-t", "--threshold", type=float,
                        # help="Probability threshold of basins",
                        # default="0.0")
    parser.add_argument('-m', '--masks', type=argparse.FileType('r'), nargs='+',
                        help="Path to files containing the state masks.",
                        required=True)
    parser.add_argument('-c', '--column', type=int, default='-1', dest='col',
                        help="Column that contains the free energy starting at 0. \
                        Defaults to -1 (the last).")
    parser.add_argument('-np', '--numprocs', type=int, default='1',
                        help="Number of parallel processes")
    parser.add_argument('-o', '--outfile',
                        help="Name of the output file(s). Default is 'delta_F'")
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


def calculate_delta_F(filename, kT, masks, col):
    """Calculates the free energy difference between states

    If more than two are specified, this returns the difference to the first state for all others

    Arguments
    ---------
    filename : path to the fes file
    kT       : energy in units of kT
    masks    : a list of boolean numpy arrays resembling the states
    col      : the file column containing the fes

    Returns
    -------
    delta_F  : a list of doubles containing the free energy difference to the first state
    """

    fes = np.genfromtxt(filename).T[col]  # assumes that last column is free energy
    if len(fes) != len(masks[0]):
        raise ValueError('Masks and FES of file {} are not of the same length ({} and {})'
                         .format(filename, len(masks[0]), len(fes)))

    probabilities = np.exp(- fes / float(kT))
    state_probs = [np.sum(probabilities[m]) for m in masks]
    delta_F = [- kT * np.log(state_probs[i]/state_probs[0]) for i in range(1, len(state_probs))]
    return delta_F


def main():
    # read in cli arguments, define constants
    args = parse_args()
    fmt_times = '%10d'
    fmt_error = '%14.9f'

    masks = []
    for m in args.masks:
        try:
            mask = np.genfromtxt(m).astype('bool')  # could also save in binary but as int/bool is more readable
        except OSError:
            print('Error: Specified masks file "{}" not found'.format(m))
            raise
        masks.append(mask)
    num_states = len(masks)
    if not all(len(m) == len(masks[0]) for m in masks):
        raise ValueError('Masks are not of the same length ({})'
                         .format(*[len(m) for m in masks]))

    if args.fd == 'f':
        print(fmt_error % calculate_delta_F(args.path, args.kT, masks, args.col))

    elif args.fd == 'd':
        folders = hlpmisc.get_subfolders(args.path)

        if not folders:  # no subdirectories found - use only given one
            if args.path[-1] != os.path.sep:
                args.path += os.path.sep  # add possibly missing "/"
            folders = [args.path]
            if args.average:
                raise ValueError("No subdirectories found. Averaging not possible.")

        # has to be done here as it's previously not clear if multiple folders are involved
        outfilenames, avgfilename = get_outfilenames(args.outfile, folders)

        files, times = hlpmisc.get_fesfiles(folders[0])

        pool = Pool(processes=args.numprocs)

        allfilenames = [os.path.join(d, f) for d in folders for f in files]
        delta_F = pool.map(partial(calculate_delta_F, kT=args.kT, masks=masks, col=args.col), allfilenames)
        delta_F = np.array(delta_F).reshape((len(folders), len(files), num_states-1))

        fields = ['time']
        if num_states == 2:
            fields.append('deltaF')
        else:
            fields += ['deltaF_1_' + str(i) for i in range(2,num_states+1)]
        header = plmdheader.PlumedHeader(fields)
        header.set_constant('kT', args.kT)
        fmt = [fmt_times] + [fmt_error]*(num_states-1)

        for i, f in enumerate(outfilenames):
            # if len(folders) == 1:  # delta_F contains only one dataset
                # outdata = np.vstack((times,delta_F.T)).T
            # else:
            outdata = np.c_[times, delta_F[i]]
            hlpmisc.backup_if_exists(f)
            np.savetxt(f, outdata, header=str(header), fmt=fmt,
                       comments='', delimiter=' ', newline='\n')

        if args.average:
            avg_delta_F = np.average(delta_F, axis=0)
            stddev = np.std(delta_F, axis=0, ddof=1)
            outdata = np.c_[times, avg_delta_F, stddev]
            fields = ['time']
            if num_states == 2:
                fields += ['deltaF.avg', 'deltaF.stddev']
            else:
                fields += ['deltaF_1_' + str(i) + '.avg' for i in range(2,num_states+1)]
                fields += ['deltaF_1_' + str(i) + '.stddev' for i in range(2,num_states+1)]
            header.fields = fields
            fmt += [fmt_error]
            hlpmisc.backup_if_exists(avgfilename)
            np.savetxt(avgfilename, outdata, header=str(header),
                       fmt=fmt, comments='', delimiter=' ', newline='\n')


if __name__ == '__main__':
    main()

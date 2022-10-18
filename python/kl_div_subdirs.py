#!/usr/bin/env python3
"""
Calculates Kullback-Leibler divergence of a set of FES files with respect to some reference
The set of files is in numbered subdirs
"""

import argparse
from functools import partial
import os
from multiprocessing import Pool
import numpy as np
from helpers import misc as hlpmisc
from kl_div import kl_div_to_ref, fes_to_prob


def parse_args():
    """Get cli args"""
    parser = argparse.ArgumentParser()
    parser.add_argument('path',
                        help="Path to the folder to be evaluated")
    parser.add_argument('-s','--subfolders', nargs='+',
                        help="Manually specify the subfolders to evaluate")
    parser.add_argument('-r', '--ref',
                        help="Path to the reference FES file",
                        required=True)
    parser.add_argument('-k', '--kT', type=float,
                        help="Value of kT for the FES files (in matching units)",
                        required=True)
    parser.add_argument('-d', '--dim', type=int, default=1,
                        help="Dimension of the FES files, defaults to 1")
    parser.add_argument('-o', '--outfile', default="kl_div",
                        help="Name of the output file(s), defaults to \"kl_div\"")
    parser.add_argument("-a", "--average", action='store_true',
                        help="Also calculate average over runs.")
    parser.add_argument('-i', '--invert', action='store_true',
                        help="Use the reference probabilities as Q(x) instead of P(x)")
    parser.add_argument("-np", "--numprocs", type=int, default="1",
                        help="Number of parallel processes")
    args = parser.parse_args()
    return args




if __name__ == '__main__':
    args = parse_args()

    # read reference fes file
    try:
        ref = np.genfromtxt(args.ref).T
    except IOError:
        print("Reference file not found!")
    if ref.shape[0] != args.dim + 1: # dim colvar columns + data colum
        raise ValueError("Specified dimension and dimension of reference file do not match.")
    ref = fes_to_prob(ref[args.dim], args.kT) # overwrite ref with the probabilities

    # get folders and files
    folders = args.subfolders
    if not folders:
        folders = hlpmisc.get_subfolders(args.path)
        if not folders:  # is empty
            print("There are no subfolders of the form '[0-9]*' at the specified path.")
            if args.average:
                raise ValueError("Averaging not possible. Are you sure about the -a option?")
            else:
                print("Using only the FES files of the base directory.")
                folders = [args.path]
    for f in folders:
        if not os.path.isdir(f):
            raise ValueError("Could not find specified subfolder {}.".format(f))

    files, times = hlpmisc.get_fesfiles(folders[0])  # assumes all folders have the same files
    if not files:
        raise ValueError("No FES files found in first subdir. Are you sure the path is correct?")

    allfilenames = [os.path.join(folder, f) for folder in folders for f in files]

    pool = Pool(processes=args.numprocs)

    kl = pool.map(partial(kl_div_to_ref, kT=args.kT, ref=ref, dim=args.dim, inv=args.invert), allfilenames)
    kl = np.array(kl).reshape(len(folders),len(files)) # put in matrix form

    fileheader =     "#! FIELDS time kl_div"
    fileheader += ("\n#! SET kT " + str(args.kT))

    for i, folder in enumerate(folders):
        outfile = os.path.join(folder, args.outfile)
        hlpmisc.backup_if_exists(outfile)
        np.savetxt(outfile, np.vstack((times, kl[i])).T, header=fileheader,
                   delimiter=' ', newline='\n')

    avgheader =     "#! FIELDS time kl_div stddev"
    avgheader += ("\n#! SET kT " + str(args.kT))

    if args.average:
        avgfile = os.path.join(args.path, args.outfile) # in base dir
        hlpmisc.backup_if_exists(avgfile)
        avg_kl = np.average(kl, axis=0)
        stddev = np.std(kl, axis=0, ddof=1)
        np.savetxt(avgfile, np.vstack((times, avg_kl, stddev)).T, header=avgheader,
                   delimiter=' ', newline='\n')

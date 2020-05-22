import argparse
import os
from multiprocessing import Pool
import numpy as np
from helpers import misc as hlpmisc


def parse_args():
    """Get cli args"""
    parser = argparse.ArgumentParser()
    parser.add_argument('path',
                        help="Path to the folder to be evaluated")
    parser.add_argument("-np", "--numprocs", type=int, default="1",
                        help="Number of parallel processes")
    args = parser.parse_args()
    return args

def get_min_pos(filename):
    """Get minimum position of FES"""
    fes, colvar = np.genfromtxt(filename).T[0:2]
    return colvar[np.where(fes == np.nanmin(fes))][0] # only first minimum...


def main():
    args = parse_args()

    fmt_times = '%10d'
    fmt_colvar= '%14.9f'

    folders = hlpmisc.get_subfolders(args.path)
    files, times = hlpmisc.get_fesfiles(folders[0]) # assumes all folders have the same files
    paths = [os.path.join(d, f) for d in folders for f in files]

    # pool = Pool(processes=args.numprocs)
    min_positions = []
    for p in paths:
        min_positions.append(get_min_pos(p))
    min_positions= np.array(min_positions).reshape(len(folders),len(files)) # put in matrix form

    fmt = [fmt_times] + [fmt_colvar] * len(folders)
    np.savetxt("min_positions", np.vstack((times,min_positions)).T, fmt=fmt)


if __name__ == '__main__':
    main()

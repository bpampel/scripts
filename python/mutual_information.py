#!/usr/bin/env python3
"""
Calculate the mutual information between 2 or more trajectories
"""

import argparse
import numpy as np
import itertools


def parse_args():
    """Get cli args"""
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+',
                        help="Path to the files to analyze")
    parser.add_argument('-c', '--column', type=int, default=2,
                        help="Column of files to analyze, defaults to 2")
    parser.add_argument('-b', '--bins', type=int, default=10,
                        help="Number of bins for histogram")
    parser.add_argument('-r', '--hist-range', dest='hist_range', type=float, nargs=2,
                        help="Number of bins for histogram")
    args = parser.parse_args()
    # exclude some options that don't make sense together

    return args


def calc_mutual_information(x, y, bins, ra, ix=None, iy=None):
    """Calculate mutual information between vectors x and y

    Arguments
    ---------
    x, y : vectors with data points
    bins : bins for histogramming the vectors
    ra   : range of the histogram
    ix   : identifier of x when printing out
    iy   : identifier of y when printing out

    Returns
    -------
    mutual_information: mutual information between x and y
    """
    c_xy = np.histogramdd(np.c_[x, y], bins, [ra, ra])[0]
    c_x = np.histogram(x,bins, ra)[0]
    c_y = np.histogram(y,bins, ra)[0]

    h_x = shannon_entropy(c_x)
    h_y = shannon_entropy(c_y)
    h_xy = shannon_entropy(c_xy)

    mutual_information = h_x + h_y - h_xy
    print(f"{ix}, {iy}: {h_x} + {h_y} - {h_xy} = {mutual_information}")
    return mutual_information, h_xy


def shannon_entropy(c):
    """Get shannon entropy of histogram data"""
    c_normalized = c / float(np.sum(c))
    c_normalized = c_normalized[np.nonzero(c_normalized)]
    h = -sum(c_normalized* np.log2(c_normalized))
    return h


def main():
    # read in cli arguments, define constants
    args = parse_args()
    fmt_times = '%10d'
    fmt_mi = '%14.9f'

    if len(args.files) < 2:
        raise ValueError("At least 2 files must be given")

    data = [np.genfromtxt(f)[1000000:,args.column-1] for f in args.files]

    n = len(data)
    mut_info = np.zeros((n,n))  # need only half the matrix but create full one
    shared_entropy = np.zeros((n,n))
    # use i1 and i2 as indices to know which elements are currently processed
    print(f"x, y: entr_x + entr_y - entr_xy = mutual_information")
    for i1, i2 in itertools.combinations(range(len(data)), 2):
        mut_info[i1, i2], shared_entropy[i1,i2] = calc_mutual_information(data[i1], data[i2], args.bins, args.hist_range, i1, i2)

    norm_mut_info = 1 - (mut_info / shared_entropy)

    print("\nMutual information:")
    print(mut_info)


    print("\nNormalized mutual information:")
    print(norm_mut_info)

if __name__ == '__main__':
    main()

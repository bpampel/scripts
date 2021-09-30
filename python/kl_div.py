#!/usr/bin/env python3
"""Script hosting the basic KL divergence math. Can be used standalone to get some simple KL div to reference"""

import argparse
import numpy as np
import sys


# first the general math functions that can be imported in custom scripts
def kl_div(p, q):
    """
    Calculates the Kullback-Leibler divergence of two probability distributions

    Arguments
    ---------
    p : numpy array containing the data / reference probabilities
    q : numpy array containing the model probabilities

    Returns
    -------
    kl_div : a double containing the KL divergence
    """
    x = np.where((p > 0) & (q > 0), p * np.log(p/q), 0)  # filter out log(0)

    return np.sum(x)


def kl_div_to_ref(filename, ref, kT, dim, inv=False):
    """
    Calculates the Kullback-Leibler divergence of a FES file to a reference

    Arguments
    ---------
    filename : path to the fes file
    ref      : numpy array holding the reference probabilities
    kT       : energy in units of kT

    Returns
    -------
    kl_div : a double containing the KL divergence
    """
    fes = np.genfromtxt(filename).T[dim]
    if fes.shape != ref.shape:
        raise ValueError("Number of elements of reference and file " + filename + " does not match")
    prob = fes_to_prob(fes, kT)

    if inv:
        return kl_div(prob, ref)
    else:
        return kl_div(ref, prob)


def fes_to_prob(fes, kT):
    """Returns normalized probability distribution from free energy surface"""
    prob = np.exp(- fes / float(kT))
    return normalize_distribution(prob)


def normalize_distribution(x):
    """Normalize a probability distribution"""
    norm = 1 / np.sum(x)
    return x * norm


# now the code for standalone execution
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('path', type=str, nargs='+', help="Path to file(s)")
    parser.add_argument('-r', '--ref', required=True,
                        type=argparse.FileType('r', encoding='UTF-8'),
                        help="Path to reference")
    parser.add_argument('-kT', '--temp', type=float, dest='kt', required=True,
                        help="Energy (in units of kT) of the FES file")
    return parser.parse_args()


def main():
    args = parse_args()
    try:
        ref = np.genfromtxt(args.ref).T
    except IOError:
        sys.exit("Reference file could not be parsed correctly!")
    dim = ref.shape[0] - 1
    colvar = ref[0:dim]
    ref = kl_div.fes_to_prob(ref[dim], args.kt)

    # calculate all
    fes_files = args.path
    kl_values = [kl_div.kl_div_to_ref(f, ref, args.kt, dim) for f in fes_files]

    # output: simply print to screen
    print('kl_div   filename')
    for f, kl_val in zip(fes_files, kl_values):
        print(f"{kl_val:.10f}   {f}")

if __name__ == '__main__':
    main()

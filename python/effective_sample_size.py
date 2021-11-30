#!/usr/bin/env python3
"""Calculate the effective sample size of a biased simulation from the weights"""

import argparse
import numpy as np


def parse_args():
    """Get cli args"""
    parser = argparse.ArgumentParser()
    parser.add_argument("path",
                        help="Path to the colvar file")
    parser.add_argument("-c", "--biascol", type=int,
                        help="Column of the file containing the bias values")
    parser.add_argument("-t", "--thermalenergy", type=float,
                        help="Thermal energy of the simulation (kT) in matching units")
    return parser.parse_args()


def calc_eff_sample_size(weights):
    """Calculate the effective sample size

    Formula from eq (13) of Invernizzi, Piaggi, Parrinello, Phys Rev X (2020)
    """
    return np.sum(weights)**2 / np.sum(weights**2)


def calc_normalized_weights(biasvals, kt):
    """Calculate normalized weights from bias values"""
    biasmax=np.max(biasvals)
    return np.exp((biasvals-biasmax)/kt)

def main():
    args = parse_args()
    biasvals = np.genfromtxt(args.path)[:,args.biascol-1] # python cols start from 0
    weights = calc_normalized_weights(biasvals, args.thermalenergy)
    eff_sample_size = calc_eff_sample_size(weights)
    print(eff_sample_size)


if __name__ == '__main__':
    main()

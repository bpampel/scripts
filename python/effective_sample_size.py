#!/usr/bin/env python3
"""Calculate the effective sample size of a biased simulation from the weights"""

import argparse
import numpy as np


def parse_args():
    """Get cli args"""
    parser = argparse.ArgumentParser()
    parser.add_argument("path", nargs='+',
                        help="Path to the colvar file(s)")
    parser.add_argument("-c", "--biascol", type=int,
                        help="Column of the file containing the bias values")
    parser.add_argument("-t", "--thermalenergy", type=float,
                        help="Thermal energy of the simulation (kT) in matching units")
    return parser.parse_args()


def eff_sample_size_from_weights(weights):
    """Calculate the effective sample size

    Formula from eq (13) of Invernizzi, Piaggi, Parrinello, Phys Rev X (2020)
    """
    return np.sum(weights) ** 2 / np.sum(weights ** 2)


def calc_normalized_weights(biasvals, kt):
    """Calculate normalized weights from bias values"""
    biasmax = np.max(biasvals)
    return np.exp((biasvals - biasmax) / kt)


def analyze_simulation(filepath, biascol, kt):
    """Calculate the effective sample size and statistics of the weights

    param filepath: path to colvar file, or list of filepaths for multiple walkers
    param biascol: column of file containing the bias values
    param kt: thermal energy of simulation
    returns: (effective_sample_size, weights_mean, weights_stddev)
    """
    if isinstance(filepath, list):
        biasvals = []
        for f in filepath:
            biasvals.append(np.genfromtxt(f)[:, biascol])
        biasvals = np.concatenate(biasvals)
    else:
        biasvals = np.genfromtxt(filepath)[:, biascol]
    weights = calc_normalized_weights(biasvals, kt)
    weights_mean = np.mean(weights)
    weights_stddev = np.std(weights, ddof=1)
    eff_sample_size = eff_sample_size_from_weights(weights)
    return (eff_sample_size, weights_mean, weights_stddev)


def main():
    args = parse_args()
    sample_size, w_mean, w_stddev = analyze_simulation(
        args.path,
        args.biascol - 1,  # python cols start from 0
        args.thermalenergy,
    )
    print(f"{sample_size}    {w_mean}    {w_stddev}")


if __name__ == "__main__":
    main()

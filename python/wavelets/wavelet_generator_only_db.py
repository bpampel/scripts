#!/usr/bin/env python3
"""
    Calculate the scaling function and its derivative for Daubechies Wavelets
    Can also be used to print out the filter coefficients

    Copyright (C) 2021 Benjamin Pampel

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import argparse
import numpy as np
from scipy.special import comb


def parse_args():
    """Get cli args"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-N", "--moments",
                        required=True, type=int,
                        help="Number of vanishing moments of the wavelet",
    )
    parser.add_argument("-r", "--recursion",
                        type=int, default=6,
                        help="Depth of recursion,\
                              i.e. 2**d points will be calculated per integer.\n\
                              Defaults to 6 (64 points per int)."
    )
    parser.add_argument("-d", "--derivs",
                        type=int, default=0,
                        help="Number of derivatives to also calculate.\n\
                              Defaults to 0 (no derivatives)."
    )
    parser.add_argument("-f", "--filename",
                        help="Name of the output file.\n\
                              Defaults to DbN.data if wavelets are calculated.",
    )
    parser.add_argument("-c", "--coeffs",
                        action="store_true",
                        help="Only calculate the filter coefficients.\n\
                              Output is by default to screen\
                              but can be to file if the -f flag is given.",
    )
    parser.add_argument("-n", "--coeffsnorm",
                        type=float, default=np.sqrt(2),
                        help="Normalization of the filter coefficients to be printed.\n\
                              Will be ignored if the -c flag is not specified.\n\
                              Defaults to sqrt(2).",
    )
    return parser.parse_args()


def filter_coeffs(N, norm=np.sqrt(2)):
    """Calculate the filter coefficients for Daubechies Wavelets

    see Strang & Nguyen - "Wavelets and Filters", 1997, ch. 5.5
    :param p: number of vanishing moments
    :param norm: specifies the sum of the final coefficients
    """
    if N < 1:
        raise ValueError("The Wavelets must have at least 1 vanishing moment")
    if N > 34:
        raise ValueError("Sorry, the implementation does not work for N > 34")
    if norm <= 0:
        raise ValueError("The norm must be larger than 0")
    poly = np.polynomial.polynomial

    # find roots of the defining polynomial B(y)
    B_y = [comb(N + i - 1, i, exact=True) for i in range(N)]
    y = poly.polyroots(B_y)
    # calculate the roots of C(z) from the roots y via a quadratic formula
    roots = [poly.polyroots([1, 4 * yi - 2, 1]) for yi in y]
    # take the ones inside the unit circle and add roots at -1
    z = [root for pair in roots for root in pair if np.abs(root) < 1]
    z += [-1] * N
    # put together the polynomial C(z) and normalize the coefficients
    C_z = poly.polyfromroots(z)
    C_z = np.real(C_z)  # imaginary part may be non-zero because of rounding errors
    C_z *= norm / sum(C_z)
    return C_z[::-1]


def highpass_from_lowpass(h):
    """Returns the lowpass filter coefficients g from the highpass filter
    coefficients h"""
    alternating_sign_list = [(1 if i % 2 == 0 else -1) for i, _ in enumerate(h)]
    g = h[::-1] * alternating_sign_list
    return g


def m_matrices(h):
    """Sets up and returns the matrices M0 and M1 (in a list)

    The sum over the filter coefficients h should equal 1.
    Works also with the lowpass filter to calculate the wavelet.

    :param h: filter coefficients, with norm 1
    """
    # initialize matrices to zero
    n = len(h) - 1
    M0 = np.zeros((n, n), dtype=np.float64)
    M1 = np.zeros((n, n), dtype=np.float64)
    # 'c style' construction
    for i in range(n):
        for j in range(n):
            if 0 <= 2 * i - j <= n:
                M0[i, j] = 2 * h[2 * i - j]
            if -1 <= 2 * i - j <= n - 1:
                M1[i, j] = 2 * h[2 * i - j + 1]
    return [M0, M1]


def normalize_eigenvector(eigvec, derivnum):
    """Normalize the eigenvector

    Required to obtain the values of phi (or its derivatives) at the integer values.
    :param eigvec: eigenvector of the M0 matrix
    :param derivnum: derivative currently calculated
    :return: phi at the integer values
    """
    weighted_sum = np.sum(eigvec * ((-np.arange(len(eigvec))) ** derivnum))
    norm = np.math.factorial(derivnum) / weighted_sum
    return eigvec * norm


def wavelet(N, d=6, derivs=0):
    """Calculate the scaling and wavelet function for Daubechies Wavelets

    :param N: number of vanishing moments N
    :param d: recursion number. returned array will have 2**d values per integer
    :param derivs: number of derivatives to also calculate
    :return (x, phi, psi): phi is the scaling and psi the wavelet function
                           both are lists containing the function and requested derivatives
    """
    if derivs > N - 1:
        raise ValueError(f"Only {N-1} derivatives exist but requested {derivs}")

    h = filter_coeffs(N, 1)
    g = highpass_from_lowpass(h)

    H = m_matrices(h)
    G = m_matrices(g)

    step = 1 << d  # number of values between integers
    # set up arrays of values to be calculated
    phi = np.empty((N, (2 * N - 1) * step), dtype=np.float64)
    psi = np.empty((N, (2 * N - 1) * step), dtype=np.float64)

    # get eigenvalues and vectors of matrix
    H0_eigvals, H0_eigvecs = np.linalg.eig(H[0])

    # identify the indices corresponding to the required derivatives in H0_eigvals and H0_eigvecs
    # the eigenvalues are 2**(-deriv) starting with 0
    # prefix 'dy' is for the dyadic values
    dy_eigval_indices = [
        np.argmin(np.absolute(H0_eigvals - 2 ** (-j))) for j in range(derivs + 1)
    ]

    values_at_int = np.empty((N, 2 * N - 1), dtype=np.float64)

    for j, k in enumerate(dy_eigval_indices):
        # j is order of derivative (0:derivs), k the position of the corresponding eigenvector
        values_at_int[j] = normalize_eigenvector(H0_eigvecs[:, k], j)

        # multiply matrices with factor (less flops)
        factor = 1 << j
        H_temp = [factor * H[0], factor * H[1]]
        G_temp = [factor * G[0], factor * G[1]]

        # fill first two datasets by hand
        binarydict = {"0": values_at_int[j]}
        binarydict["1"] = H_temp[1] @ values_at_int[j]
        phi[j][::step] = binarydict["0"]
        phi[j][step >> 1 :: step] = binarydict["1"]
        psi[j][::step] = G_temp[0] @ values_at_int[j]
        psi[j][step >> 1 :: step] = G_temp[1] @ values_at_int[j]

        # do the recursion
        oldbits = ["1"]
        for depth in range(2, d + 1):
            newbits = ["%d%s" % (new, old) for new in [0, 1] for old in oldbits]
            for binary in newbits:
                start = int(binary, 2) * step >> depth
                firstbit = int(binary[0])
                binarydict[binary] = H_temp[firstbit] @ binarydict[binary[1:]]
                phi[j][start::step] = binarydict[binary]
                psi[j][start::step] = G_temp[firstbit] @ binarydict[binary[1:]]
            oldbits = newbits

    # corresponding x values
    x = np.arange(0, 2 * N - 1, 1 / step)

    return (x, phi, psi)


if __name__ == "__main__":
    args = parse_args()

    if args.coeffs:  # print only file coeffs
        coeffs = filter_coeffs(args.moments, args.coeffsnorm)
        if args.filename is None:  # print to screen
            [print("{:.32e}".format(c)) for c in coeffs][0]
        else:
            header = "#! FIELDS coeff"
            header += "\n#! SET type Db" + str(args.moments)
            np.savetxt(args.filename, coeffs.T, header=header)

    else:  # do the main program, i.e. calculate and print the wavelet functions
        x, scaling_func, wavelet = wavelet(args.moments, args.recursion, args.derivs)
        if args.filename is None:  # default names
            args.filename = f"Db{args.moments}.data"
        header = "#! FIELDS x Phi Psi"
        header += "".join([f" Phi_d{i} Psi_d{i}" for i in range(1, args.derivs + 1)])
        header += "\n#! SET type Db" + str(args.moments)

        # collect all data into one array
        data = np.empty((3 + args.derivs * 2, len(x)))
        data[0] = x
        for deriv in range(args.derivs + 1):
            data[2 * deriv + 1] = scaling_func[deriv]
            data[2 * deriv + 2] = wavelet[deriv]
        np.savetxt(args.filename, data.T, header=header)

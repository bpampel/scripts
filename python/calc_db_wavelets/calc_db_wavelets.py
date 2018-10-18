#!/usr/bin/env python3
'''
Calculate the scaling function and its derivative for Daubechies Wavelets
'''

import sys
from numpy.polynomial import polynomial as poly
import numpy as np
from scipy.special import comb

def db_filter_coeffs(p, norm=1):
    ''' Calculate the filter coefficients for Daubechies Wavelets

    c.f. Strang, Nguyen - "Wavelets and Filters" ch. 5.5
    :param p: number of vanishing moments
    :param norm: specifies the euclidean norm of the final coefficients
    '''
    if p < 1:
        raise ValueError("The order of the Wavelets must be at least 1")
    if p > 34:
        raise ValueError("Sorry, the method does currently not work for orders higher than 34")
    if norm <= 0:
        raise ValueError("The norm must be larger than 0")
    # find roots of the defining polynomial B(y)
    B_y = [comb(p + i - 1, i, exact=1) for i in range(p)]
    y = poly.polyroots(B_y)
    # calculate the roots of C(z) from the roots y via a quadratic formula
    z = [[root for root in poly.polyroots([1, 4*yi - 2, 1]) if abs(root) < 1][0] for yi in y]
    z += [-1]*p
    # put together the polynomial C(z) and normalize the coefficients
    C_z = np.real(poly.polyfromroots(z))
    C_z *= norm * np.sqrt(2) / sum(C_z)
    return C_z[::-1]


def highpass_from_lowpass(h):
    '''Returns the lowpass filter coefficients g from the highpass filter
       coefficients h'''
    alternating_sign_list = [(1 if i%2 == 0 else -1) for i, _ in enumerate(h)]
    g = h[::-1] * alternating_sign_list
    return g


def m_matrices(h):
    '''Sets up and returns the matrices M0 and M1 (in a list)
       that are needed to get the values of the scaling function
       at integer points.
       The sum over the filter coefficients h should equal 1.
       Works also with the lowpass filter to calculate the wavelet.
    '''
    # initialize matrices to zero
    n = len(h)-1
    M0 = np.zeros((n, n), dtype=np.float64)
    M1 = np.zeros((n, n), dtype=np.float64)
    # 'c style' but better for transfering anyway
    for i in range(n):
        for j in range(n):
            if 0 <= 2*i-j <= n:
                M0[i,j] = 2*h[2*i -j]
            if -1 <= 2*i-j <= n - 1:
                M1[i,j] = 2*h[2*i -j+1]
    # print("M0:\n{}\n".format(M0))
    # print("M1:\n{}\n".format(m1/2))
    return [M0, M1]


def normalize_eigenvector(eigvec, derivnum):
    '''Normalize the eigenvector to retrieve the values of phi (or its
       derivatives) at the integer values.'''
    weighted_sum = np.sum(eigvec * ((-np.arange(len(eigvec))) ** derivnum))
    # print(weighted_sum)
    norm = np.math.factorial(derivnum) / weighted_sum
    # print(norm)
    return eigvec * norm


def db_wavelet(p, d=6):
    '''Calculate the scaling and wavelet function for Daubechies Wavelets

    :param p: order p of (equivalent to the number of vanishing moments)
    :param d: recursion number. returned array will have 2**d values per integer
    :return: (x, phi, psi) where phi is the scaling and psi the wavelet function
    '''
    n = 2*p - 1
    h = calc_filter_coeffs(p, 1/np.sqrt(2))
    g = calc_highpass_from_lowpass(h)
    # print("h: {}\n".format(h))
    H = setup_m_matrices(h)
    G = setup_m_matrices(g)

    step = 1 << d # number of values between integers
    # set up arrays of values to be calculated
    phi = np.empty((p, (2*p-1) * step), dtype=np.float64)
    psi = np.empty((p, (2*p-1) * step), dtype=np.float64)

    # get eigenvalues of matrices
    H0_eigvals, H0_eigvecs = np.linalg.eig(H[0])
    # This could be simplified: We already know the eigenvalues
    # --> simply find the corresponding vectors

    # prefix 'dy' is for the dyadic values
    dy_eigval_indices = [np.argmin(np.absolute(H0_eigvals - 2**(-j))) for j in range(p)]
    values_at_int = np.empty((p,n), dtype=np.float64)
    for j, k in enumerate(dy_eigval_indices):
        # j is order of derivative (0:p-1), k the position of the corresponding eigenvector
        values_at_int[j] = normalize_eigenvector(H0_eigvecs[:,k], j)

        # print("Deriv {}: (factorial check: {})\n{}\n".format(
            # j, np.sum((-np.arange(n))**j* values_at_int[j]),
            # values_at_int[j])
             # )

        # multiply matrices with factor (less flops)
        factor = 1 << j
        H_temp = [factor * H[0], factor * H[1]]
        G_temp = [factor * G[0], factor * G[1]]

        # fill first two datasets by hand
        binarydict = {'0': values_at_int[j]}
        binarydict['1'] = H_temp[1] @ values_at_int[j]
        phi[j][::step] = binarydict['0']
        phi[j][step >> 1::step] = binarydict['1']
        psi[j][::step] = G_temp[0] @ values_at_int[j]
        psi[j][step >> 1::step] = G_temp[1] @ values_at_int[j]
        if j == 0:
            print("Values at ints: \n{}\n".format(G_temp[0] @ values_at_int[j]))

        oldbits = ['1']
        for depth in range(2, d + 1):
            # is there a nicer way?
            newbits = ['%d%s' % (new, old) for new in [0, 1] for old in oldbits]
            for binary in newbits:
                start = int(binary, 2) * step>>depth
                firstbit = int(binary[0])
                binarydict[binary] = H_temp[firstbit] @ binarydict[binary[1:]]
                phi[j][start::step] = binarydict[binary]
                psi[j][start::step] = G_temp[firstbit] @ binarydict[binary[1:]]
            oldbits = newbits

    # corresponding x values
    x = np.arange(0, 2*p - 1, 1 / step)

    return (np.vstack((x, phi)), np.vstack((x, psi)))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        order = int(input("Please specify the order of the Daubechies Wavelets to "
                          "be calculated (1-34):\n"))
    else:
        order = int(sys.argv[1])

    scaling_func, wavelet = calc_scaling_function(order, 6)
    np.savetxt("Phi_{}.dat".format(order), np.transpose(scaling_func))
    np.savetxt("Psi_{}.dat".format(order), np.transpose(wavelet))

#!/usr/bin/env python3
'''
Calculate the filter coefficients of Daubechies Wavelets
'''

import sys
from numpy.polynomial import polynomial as poly
import numpy as np
from scipy.special import comb


def calc_filter_coeffs(p, norm=1):
    '''Calculate the filter coefficients for Daubechies Wavelets
       c.f. Strang, Nguyen - "Wavelets and Filters" ch. 5.5
       p is the number of vanishing moments
       norm specifies the euclidean norm of the final coefficients'''
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


def setup_m_matrices(h):
    '''Sets up and returns the matrices m0 and m1 (in a list)
    that are needed to get the values of the scaling function
    at integer points'''
    # initialize matrices to zero
    n = len(h)-1
    m0 = np.zeros((n, n), dtype=np.float64)
    m1 = np.zeros((n, n), dtype=np.float64)
    # 'c style' but better for transfering anyway
    for i in range(n):
        for j in range(n):
            if 0 <= 2*i-j <= n:
                m0[i,j] = 2*h[2*i -j]
            if -1 <= 2*i-j <= n - 1:
                m1[i,j] = 2*h[2*i -j+1]
    print("m0:\n{}\n".format(m0/2))
    print("m1:\n{}\n".format(m1/2))
    return [m0, m1]


def calc_scaling_function(p, d=8):
    '''Calculate the scaling function Phi(x) of Daubiches Wavelets
       for a given order p with 2^d values between each integer
    '''
    h = calc_filter_coeffs(p, 1/np.sqrt(2))
    print("h: {}\n".format(h))
    m = setup_m_matrices(h)
    # get eigenvalues of m0 and round them (because of precision errors)
    m0_eigvals, m0_eigvecs = np.linalg.eig(m[0])
    m0_eigvals = [round(x,8) for x in m0_eigvals]
    intvalues = m0_eigvecs[:,m0_eigvals.index(1)]
    # This could be simplified: We already know that 1 is an eigenvalue
    # --> simply find the corresponding vector

    intsum = np.sum(intvalues)
    if intsum < 0:
        intvalues *= -1
        intsum *= -1
    intvalues /= intsum

    print("Integer_phi: {}".format(intvalues))
    print("Multiplication test: {}".format((m[0] @ intvalues) - intvalues))

    step = 1 << d
    # set up array of values to be calculated
    phi = np.zeros((2*p-1) * step, dtype=np.float64)

    # fill first two datasets by hand
    binarydict = {'0': intvalues}
    binarydict['1'] = m[1] @ intvalues
    phi[::step] = binarydict['0']
    phi[step >> 1::step] = binarydict['1']

    oldbits = ['1']
    for depth in range(2, d + 1):
        # is there a nicer way?
        newbits = ['%d%s' % (new, old) for new in [0, 1] for old in oldbits]
        for binary in newbits:
            start = int(binary, 2) * step>>depth
            firstbit = int(binary[0])
            binarydict[binary] = m[firstbit] @ binarydict[binary[1:]]
            phi[start::step] = binarydict[binary]
        oldbits = newbits
    x = np.arange(0, 2*p - 1, 1 / step)

    return x, phi


if __name__ == '__main__':
    if len(sys.argv) != 2:
        order = int(input("Please specify the order of the Daubechies Wavelets to "
                      "be calculated (1-10):\n"))
    else:
        order = int(sys.argv[1])

    scaling_func = calc_scaling_function(order)
    # test it by comparison with the scipy implementation
    # from scipy.signal import wavelets as wv
    # scipyref = wv.cascade(wv.daub(order), 8)
    # differences = (scaling_func[0] - scipyref[0], scaling_func[1] - scipyref[1])
    # print("Max Delta_x: {}\nMax Delta_Phi: {}\n".format(max(differences[0]), max(differences[1])))
    np.savetxt("Phi_{}.dat".format(order), np.transpose(np.asmatrix(scaling_func)))


else:
    exit(-1)

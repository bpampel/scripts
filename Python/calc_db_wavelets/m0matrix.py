#!/usr/bin/env python3
'''
Calculate the scaling function of Daubechies Wavelets
'''


import sys
from math import sqrt
import numpy as np


def setup_h(order):
    '''Get coefficients from file'''
    return np.genfromtxt("db"+str(order)+"_coeffs")


if __name__ == '__main__':

    if len(sys.argv) != 2:
        order = input("Please specify the order of the Daubechies Wavelets to "
                      "be calculated (1-10):\n")
    else:
        order = sys.argv[1]

    h = setup_h(order)/sqrt(2)
    num_int_values = len(h) - 1
    m0 = np.zeros((num_int_values, num_int_values), dtype=np.float64)
    m1 = np.zeros((num_int_values, num_int_values), dtype=np.float64)

    # 'c style' but better for transfering anyway
    for i in range(num_int_values):
        for j in range(num_int_values):
            if 0 <= 2*i-j <= num_int_values:
                m0[i,j] = h[2*i -j]
            if -1 <= 2*i-j <= num_int_values - 1:
                m1[i,j] = h[2*i -j+1]
    m0 *= 2
    m1 *= 2
    eigvals, eigvecs = np.linalg.eig(m0)
    # np.set_printoptions(precision=3)
    # print('\nw:\n{}'.format(w))
    # print('\nv:\n{}'.format(v))

    # round eigenvalues because sometimes it is not exactly 1
    eigvals = [round(x,8) for x in eigvals]

    # xvalues = np.array(eigvecs[:,eigvals.index(1)], dtype=np.float64)
    print(eigvals)
    # mat = [m0,m1]

    # num_iterations = 10
    # for i in range(1, num_iterations):
        # temp_values = np.zeros(num_int_values * 2**i)
        # temp_values[::2] = xvalues
        # temp_values[1::2] = np.concatenate([mat[arrnum%2] @ subarray
                                            # for arrnum, subarray in
                                            # enumerate(np.split(xvalues, 2**(i-1)))])
        # xvalues = temp_values
    # np.savetxt('test.data', np.transpose(np.asmatrix(xvalues)))







    # print(h)



else:
    exit(-1)

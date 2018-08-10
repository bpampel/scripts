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
    try:
        p=int(p)
    except ValueError:
        print("Error: A non-integer value was specified as order!")
        exit(1)
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


if __name__ == '__main__':
    if len(sys.argv) != 2:
        order = input("Please specify the order of the Daubechies Wavelets to "
                      "be calculated (1-10):\n")
    else:
        order = sys.argv[1]

    print(calc_filter_coeffs(order, np.sqrt(1)))


else:
    exit(-1)

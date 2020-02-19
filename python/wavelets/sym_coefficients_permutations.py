#!/usr/bin/env python3
'''
Get all permutations of the Daubechies coefficients (to compare with the Symlet ones in literature)
Note: this includes duplicates (choosing the opposite ones gives the same coefficients reversed)
      also the roots come in quadruples, choosing the complex conjugate yields iirc also the same
'''

import sys
from numpy.polynomial import polynomial as poly
import numpy as np
from scipy.special import comb
from itertools import combinations


p=8
# find roots of the defining polynomial B(y)
B_y = [comb(p + i - 1, i, exact=True) for i in range(p)]
y = poly.polyroots(B_y)
# calculate the roots of C(z) from the roots y via a quadratic formula
roots = np.array([poly.polyroots([1, 4*yi - 2, 1]) for yi in y])

allcomb = [x for y in [[[e for e in ele] for ele in elem] for elem in [list(combinations(range(p-1),i)) for i in range(p-1)]] for x in y]
for i, comb in enumerate(allcomb):
    mask=np.array([[i in comb, i not in comb] for i in range(p-1)])
    z = roots[mask].tolist()
    z += [-1]*p
    # put together the polynomial C(z) and normalize the coefficients
    C_z = np.real(poly.polyfromroots(z)) # imaginary part may be non-zero because of rounding errors
    C_z *= (2 / sum(C_z))
    print("\n{}: {}".format(i, comb))
    [print("{:.16e}".format(i)) for i in C_z[::-1]]
    input("")


# result: num 38 is the "correct" one: choosing [0,3,4] inside

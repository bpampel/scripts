#!/usr/bin/env python2
import numpy as np
from itertools import ifilterfalse

fname='coeffs.data'
fout='norm_coeffs.data'

# for skipping headers
def iscomment(s):
    return s.startswith('#')

# filter lines that are empty
with open(fname) as f_in:
        data = filter(None, (line.rstrip() for line in f_in))

fout=open(fout,'w')

values = [] # contains values
indices = set() # contains indices

for line in ifilterfalse(iscomment,data):
    line=line.strip()
    columns=line.split()
    indices.add(int(columns[0]))
    values.append(float(columns[1]))

numi=len(indices)
valuelist = [values[x:x+numi] for x in range(0, len(values), numi)]

# now switch to numpy because it's shorter
coeffs = np.array(valuelist)
norm = np.array(valuelist[-1])

# watch out for zeros
norm[norm == 0] = 1

normcoeffs = (coeffs.T / norm[:,None]).T

print normcoeffs

np.savetxt('normcoeffs.data', np.asmatrix(normcoeffs), fmt='%1.16f', delimiter=' ', newline='\n')

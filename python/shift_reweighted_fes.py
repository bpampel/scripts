#!/usr/bin/env python2
import numpy as np
from itertools import ifilterfalse

fname='histo_ReweightDistance'
fout='histo_ReweightDistance_shifted'
temp=1 # in energy units

# for skipping headers
def iscomment(s):
    return s.startswith('#')

# filter lines that are empty
with open(fname) as fin:
        data = filter(None, (line.rstrip() for line in fin))
fin.close()

fout=open(fout,'w')

distances = []
histvalues = []

for line in ifilterfalse(iscomment,data):
    line=line.strip()
    columns=line.split()
    distances.append(float(columns[0]))
    histvalues.append(float(columns[1]))

# get minumum and shift values
logarray = -np.log(np.array(histvalues)*temp)
minvalue = min(logarray)
logarray = logarray-minvalue
print minvalue
shiftedarray = np.exp(-logarray/temp)

# now save to file, is easier in numpy (although probably inefficient)
np.savetxt(fout, np.transpose(np.matrix([np.array(distances), shiftedarray])), fmt='%1.16f', delimiter=' ', newline='\n')
fout.close()

#!/usr/bin/env python3
import numpy as np

# Total number of fes files in folder
total_files=201
# Min and max initial guesses
min_min=14
min_max=90
max_min=90
max_max=160

step=100

f = open('barrier.txt', 'w')
print('#! iter min max delta', file=f)

for i in range(total_files):
        file_name="fes.b1.iter-" + str(step*i) + ".data"
        matrix=np.genfromtxt(file_name)
        minimum=np.amin(matrix[min_min:min_max,1])
        maximum=np.amax(matrix[max_min:max_max,1])
        print(str(i), str(minimum), str(maximum), str(maximum-minimum), sep=' ', end='\n', file=f)
f.close()

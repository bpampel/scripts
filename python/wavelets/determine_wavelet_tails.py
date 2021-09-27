#!/usr/bin/env python3
'''
Calculate the scaling function and its derivative for Daubechies Wavelets
'''

import sys
import numpy as np
import wavelet_generator as wvlt

def significant_range(array, threshold):
    """Returns list with indices of first and last element above threshold"""
    above = np.nonzero(abs(array)>=threshold)
    return [above[0][0], above[0][-1]]


if __name__ == '__main__':
    for order in range(4,21):
        x, scaling_func, wavelet = wvlt.wavelet(order, 6)
        print(order,end=' ')
        for func in [scaling_func[0], wavelet[0]]:
            threshold = 0.01*np.max(func)
            indices=significant_range(func,threshold)
            pos="{} {}".format(x[indices[0]],x[indices[1]])
            print(pos,end=' ')
        print("")






        # get area above threshold
        # threshold = 0.01
        # wvlt_index = get_threshold_range(wavelet[1],threshold)
        # sclg_index = get_threshold_range(scaling_func[1],threshold)
        # print("Area with values larger than 0.01:\nWavelet [{}; {}]\nScaling function [{}; {}]".format(
                # wavelet[0][wvlt_index[0]],wavelet[0][wvlt_index[1]],scaling_func[0][sclg_index[0]],scaling_func[0][sclg_index[1]]))


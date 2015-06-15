from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
import numpy as np

data = np.array([1.,2.,3.,4.,5.])
print("data", repr(data))
bins = np.array([0., 2.5, 5.])
print("bins", repr(bins))
dig_data = np.digitize(data, bins)
h_data = np.histogram(data, bins)

data_x = np.array([1.,2.,3.,4.,5.]) # column-wise # works for histogram2d
data_y = np.array([11.,12.,13.,14.,15.]) # works for histogram2d
print("data_x", repr(data_x))
print("data_y", repr(data_y))
#data_2d = np.array([[1.,2.,3.,4.,5.], [11.,12.,13.,14.,15.]]) # column-wise # doesn't work!
data_2d = np.array([[1.,11.], [2.,12.], [3.,13.], [4.,14.], [5.,15.]]) # row-wise # works for histogramdd
print("data_2d", repr(data_2d))
#bins_2d = np.array([[0., 2.5, 5.], [10., 12.5, 15.]])
bins_2d = np.array([[0., 2., 4., 6.], [10., 12., 14., 16.]])
print("bins_2d", repr(bins_2d))
#np.digitize(data_2d, bins_2d) # doesn't work!
#np.histogram(data_2d, bins_2d) # doesn't work!
#np.histogram2d(data_2d, bins_2d) # doesn't work!
h2d_data = np.histogram2d(data_x, data_y, bins_2d)
hdd_data = np.histogramdd(data_2d, bins_2d)

import IPython; IPython.embed()
dig_data
h_data
h2d_data
hdd_data

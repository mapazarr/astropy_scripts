from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
import numpy as np
from astropy.coordinates import Angle
from gammapy.obs import ObservationGroups, ObservationGroupAxis

alt = Angle([0, 30, 60, 90], 'degree')
az = Angle([-90, 90, 270], 'degree')
ntels = np.array([3, 4])

alt_obs_group_axis = ObservationGroupAxis('ALT', alt, 'bin_edges')
az_obs_group_axis = ObservationGroupAxis('AZ', az, 'bin_edges')
ntels_obs_group_axis = ObservationGroupAxis('N_TELS', ntels, 'bin_values')

print("ALT n_bins", alt_obs_group_axis.n_bins)
print("ALT bin 2", alt_obs_group_axis.get_bin(2))
##print("ALT bins", alt_obs_group_axis.get_bins)

list_obs_group_axis = [alt_obs_group_axis, az_obs_group_axis, ntels_obs_group_axis]
array_obs_group_axis = np.array([alt_obs_group_axis, az_obs_group_axis, ntels_obs_group_axis])

obs_group = ObservationGroups(list_obs_group_axis)




##import IPython; IPython.embed()


# test the definition of observation groups

#array_ob_bins = np.array([])
#for i_axis in list_obs_group_axis:
#    array_ob_bins.append(np.arange(len(list_obs_group_axis[i_axis])))

#list_of_bins = []
#for i_axis in np.arange(len(list_obs_group_axis)):
#    list_of_bins.append(np.arange(len(list_obs_group_axis[i_axis].bins)))

array_of_bins = np.array([])
for i_axis in np.arange(len(list_obs_group_axis)):
    array_of_bins = np.append(array_of_bins, np.arange(len(list_obs_group_axis[i_axis].bins)))

##import IPython; IPython.embed()


stuff = [[0, 30, 60, 90], [-90, 90, 270], [3, 4]]
stuff = np.array([[0, 30, 60, 90], [-90, 90, 270], [3, 4]])
stuff = np.array([np.array([0, 30, 60, 90]), np.array([-90, 90, 270]), np.array([3, 4])])

grid = np.meshgrid(stuff, indexing='ij')

grid = np.meshgrid(np.array([0, 30, 60, 90]), np.array([-90, 90, 270]), np.array([3, 4]), indexing='ij')

len(grid)

##import IPython; IPython.embed()




a = np.array([11, 12])
b = np.array([21, 22])
a
b
a_ab, b_ab = np.meshgrid(a, b, indexing='ij')
a_ab
b_ab

c = np.array([31, 32])


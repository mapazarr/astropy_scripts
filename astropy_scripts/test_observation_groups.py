from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
import numpy as np
from astropy.coordinates import Angle
from astropy.table import Table
from gammapy.obs import (ObservationTable, ObservationGroups,
                         ObservationGroupAxis)
from gammapy import datasets

CACHE = 1 # set to 0 to ignore downloaded files in the cache
          # (or use rm ~/.astropy/cache)

alt = Angle([0, 30, 60, 90], 'degree')
az = Angle([-90, 90, 270], 'degree')
ntels = np.array([3, 4])

alt_obs_group_axis = ObservationGroupAxis('ALT', alt, 'bin_edges')
az_obs_group_axis = ObservationGroupAxis('AZ', az, 'bin_edges')
ntels_obs_group_axis = ObservationGroupAxis('N_TELS', ntels, 'bin_values')

print("ALT n_bins", alt_obs_group_axis.n_bins)
print("ALT bin 2", alt_obs_group_axis.get_bin(2))
print("ALT bins", alt_obs_group_axis.get_bins)
print("ALT col:", alt_obs_group_axis.to_column().meta)
print(alt_obs_group_axis.to_column())

alt_obs_group_axis_from_col = ObservationGroupAxis.from_column(alt_obs_group_axis.to_column())
print("ALT from col col:", alt_obs_group_axis_from_col.to_column().meta)
print(alt_obs_group_axis_from_col.to_column())
print("ALT from col bins", alt_obs_group_axis.bins)

print()
alt_obs_group_axis.print() # not working?
az_obs_group_axis.print() # not working?
ntels_obs_group_axis.print() # not working?
print(alt_obs_group_axis.info)
print(az_obs_group_axis.info)
print(ntels_obs_group_axis.info)

list_obs_group_axis = [alt_obs_group_axis, az_obs_group_axis, ntels_obs_group_axis]
array_obs_group_axis = np.array([alt_obs_group_axis, az_obs_group_axis, ntels_obs_group_axis])

obs_group = ObservationGroups(list_obs_group_axis)

obs_group.print_groups() # works?
obs_group.print_axes() # not working?
print(obs_group.info)
print(obs_group.obs_groups_table)

# write
obs_group_1 = obs_group
obs_group_1.write('obs_groups.ecsv')

# read
obs_group_2 = ObservationGroups.read('obs_groups.ecsv')

import IPython; IPython.embed()
assert (obs_group_1.obs_groups_table == obs_group_2.obs_groups_table).all()

print(obs_group.obs_groups_table)
print(obs_group.info)


# test columns of different lengths:
# TODO: this is an attempt for a function to pass observation axis definitions from/to a text (ECSV) file.
# It is not working because of different lengths of the arrays; possibilities:
# 1) define to_row/from_row in ObservationGroupAxis, an pack the info in a table of columns: name, unit, format, bins (the latter should be an array of either bin values or bin edges: np.array(bins))
# 2) Play with masked tables: define all columns with the same number of bins, then mask the ones that do not proceed.
# Possibility 1 can be stored in a regular CSV file, maybe?

#axis_table = Table(masked=True)
#for i_axis in list_obs_group_axis:
#    #axis_table[i_axis.name] = i_axis.bins
#    axis_table.add_column(i_axis.to_column())
#    axis_table.meta[i_axis.name + "_axis_format"] = i_axis.format
#    axis_table.meta[i_axis.name + "_axis_format"] = i_axis.format
#    print(axis_table)
#    print(axis_table.meta)
#print()
#print(axis_table.meta)
#print(axis_table)

# test apply observation groups to a dummy observation table
# using file in gammapy-extra (I also could create a dummy table)
infile = datasets.get_path('../test_datasets/obs/test_observation_table.fits',
                             location='remote', cache=CACHE)
obs_table = ObservationTable.read(infile)

print()
print(obs_table)

# wrap azimuth angles to [-90, 270) deg
# to match definition of azimuth grouping axis
obs_table['AZ'] = Angle(obs_table['AZ']).wrap_at(Angle(270., 'degree'))

obs_table_grouped = obs_group.group_observation_table(obs_table)

# wrap azimuth angles back to [0, 360) deg
obs_table['AZ'] = Angle(obs_table['AZ']).wrap_at(Angle(360., 'degree'))

print()
print(obs_table_grouped)

assert len(obs_table) == len(obs_table_grouped)
assert ((0 <= obs_table_grouped['GROUP_ID']) &
        (obs_table_grouped['GROUP_ID'] < obs_group.n_groups)).all()

print('Done.')

# TODO in ObservationGroups class:
    # add high level doc with printouts of the observation groups -> do it in the "future" inline command tool for obs groups!!!

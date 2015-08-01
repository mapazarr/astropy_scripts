from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
import numpy as np
#from astropy.table import Table, Column
from astropy.table import Table
from astropy.units import Quantity
from astropy.coordinates import Angle
from gammapy.datasets.make import make_test_observation_table
from gammapy.obs import ObservationTable

observatory_name='HESS'
n_obs = 10
#debug = True
debug = False

n_obs_start = 1

# playing with astropy table
astro_table = Table()
#col_obs_id = Column(name='obs_id', data=np.arange(n_obs_start, n_obs_start + n_obs))
#astro_table.add_column(col_obs_id)
obs_id = np.arange(n_obs_start, n_obs_start + n_obs)
astro_table['OBS_ID'] = obs_id
time_observation = Quantity(30. * np.ones_like(obs_id), 'minute')
astro_table['TIME_OBSERVATION'] = time_observation
time_live = Quantity(25. * np.ones_like(obs_id), 'minute')
astro_table['TIME_LIVE'] = time_live
print("astro_table")
print(astro_table)

print()

# get data from the table in an array
astro_array = np.array(astro_table)
print("astro_array")
print(astro_array)

print()

# convert structured (a.k.a. record) array to regular numpy array
# ref: http://stackoverflow.com/questions/5957380/convert-structured-array-to-regular-numpy-array
#astro_array = astro_array.view((float, len(astro_array.dtype.names))) # faulty (obs_id is not a float)!

# get data from the table in an array (selecting columns)
astro_array_sel = np.vstack([astro_table['OBS_ID'], astro_table['TIME_OBSERVATION'], astro_table['TIME_LIVE']]).T

#####################################################################

# playing with observation table generator
obs_table = make_test_observation_table(observatory_name, n_obs, debug=debug)

print("obs_table.summary")
print(obs_table.summary())

print()

print("obs_table.meta (header)")
print(obs_table.meta)

print()

print("obs_table")
print(obs_table)

print()

# test: assert if the length of the table is n_obs:
print("obs_table length")
print(len(obs_table))
assert len(obs_table) == n_obs

# test: assert if the TIME_START > 0:
assert (obs_table['TIME_START'] > 0).all()

# test: assert if TIME_STOP > TIME_START:
assert (obs_table['TIME_STOP'] > obs_table['TIME_START']).all()

# test: assert if RA is in the interval (0, 360) deg:
ra_min = Angle(0, 'degree')
ra_max = Angle(360, 'degree')
assert (ra_min < obs_table['RA']).all()
assert (obs_table['RA'] < ra_max).all()

# test: assert if dec is inthe interval (-90, 90) deg:
dec_min = Angle(-90, 'degree')
dec_max = Angle(90, 'degree')
assert (dec_min < obs_table['DEC']).all()
assert (obs_table['DEC'] < dec_max).all()

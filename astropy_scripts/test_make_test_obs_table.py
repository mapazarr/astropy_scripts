from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
import numpy as np
from astropy.table import Table, Column
from gammapy.datasets.make import make_test_observation_table
from gammapy.obs import ObservationTable

observatory_name='HESS'
n_obs = 10
#debug = True
debug = False

n_obs_start = 1

astro_table = Table()
col_obs_id = Column(name='obs_id', data=np.arange(n_obs_start, n_obs_start + n_obs))
astro_table.add_column(col_obs_id)
print("astro_table")
print(astro_table)

print()

obs_table = make_test_observation_table(observatory_name, n_obs, debug)

print("obs_table.info")
print(obs_table.info())

print()

print("obs_table.meta (header)")
print(obs_table.meta)

print()

print("obs_table")
print(obs_table)

print()

#test: assert if the length of the table is n_obs:
print("obs_table length")
print(len(obs_table))
assert len(obs_table) == n_obs

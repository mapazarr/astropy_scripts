import numpy as np
from astropy.table import Table, Column
from gammapy.datasets.make import generate_observation_table
from gammapy.obs import ObservationTable

observatory='HESS'
n_obs = 10

n_obs_start = 1

astro_table = Table()
col_obs_id = Column(name='obs_id', data=np.arange(n_obs_start, n_obs_start + n_obs))
astro_table.add_column(col_obs_id)
print "astro_table"
print(astro_table)

print ""

obs_table = generate_observation_table(observatory, n_obs)

print "obs_table.info"
print obs_table.info()

print ""

print "obs_table"
print(obs_table)

print ""

#test: assert if the length of the table is n_obs:
print "obs_table length"
print(len(obs_table))
assert len(obs_table) == n_obs

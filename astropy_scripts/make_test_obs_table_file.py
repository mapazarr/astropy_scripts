from gammapy.datasets import make_test_observation_table
from gammapy.obs import ObservationTable

SAVE = 0

outfile = 'test_observation_table.fits'
overwrite = True

observatory_name='HESS'
#n_obs=10
n_obs=100

observation_table = make_test_observation_table(observatory_name=observatory_name,
                                                n_obs=n_obs)

import IPython; IPython.embed() # and save manually

# write
if SAVE:
    print("Writing {}.".format(outfile))
    observation_table.write(outfile, overwrite=overwrite)

# read:
#observation_table = ObservationTable.read(outfile)

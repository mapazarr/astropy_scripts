from gammapy.datasets import make_test_observation_table
from gammapy.obs import ObservationTable

SAVE = 0

outfile = 'test_observation_table.fits'
overwrite = True

#observation_table = make_test_observation_table('HESS', 10)
observation_table = make_test_observation_table('HESS', 100)

import IPython; IPython.embed() # and save manually

# write
if SAVE:
    print("Writing {}.".format(outfile))
    observation_table.write(outfile, overwrite=overwrite)

# read:
#observation_table = ObservationTable.read(outfile)

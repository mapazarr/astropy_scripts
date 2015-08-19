from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
import os
import numpy as np
from gammapy.datasets import make_test_dataset

CLEAN_WORKING_DIR = 1 # remove existing dataset dir

observatory_name = 'HESS'
#n_obs = 10
n_obs = 2
datestart = None
dateend = None
#random_state = 'random-seed'
#random_state = 0 # this is not deterministic!!!
random_state = np.random.RandomState(seed=0)

fits_path = '/home/mapaz/astropy/development_code/astropy_scripts/astropy_scripts/' + 'test_dataset'
overwrite = False
#overwrite = True

# Need to make sure the working dir is clean, otherwise old
# files could be mixed up in the new models!
if CLEAN_WORKING_DIR:
    print("Cleaning working dir.")
    command = "rm test_dataset -fr"
    print(command)
    os.system(command)

make_test_dataset(fits_path=fits_path,
                  overwrite=overwrite,
                  observatory_name=observatory_name,
                  n_obs=n_obs,
                  datestart=datestart,
                  dateend=dateend,
                  random_state=random_state)

# test number of files created
n_event_list_files = sum(len([f for f in fs if f.lower().endswith('.fits.gz')])
                         for _, _, fs in os.walk(fits_path))
assert n_event_list_files == n_obs

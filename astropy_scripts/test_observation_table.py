from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
from astropy.coordinates import Angle
from gammapy.obs import DataStore, ObservationTable
from gammapy.datasets.make import make_test_observation_table

# case 1: using `~gammapy.obs.DataStore.make_observation_table` (works)
#         this function returns an astropy Table object
HESSFITS_MPP = '/home/mapaz/astropy/gammapy_tutorial/HESS_fits_data/pa/Model_Deconvoluted_Prod26/Mpp_Std/'
data_store = DataStore(dir=HESSFITS_MPP)
obs_table1 = data_store.make_observation_table()
#observatory_name='HESS'
#n_obs = 100
#observation_table = make_test_observation_table(observatory_name, n_obs)
zenith_min = Angle(20., 'degree')
zenith_max = Angle(40., 'degree')
zenith = Angle(90., 'degree') - obs_table1['ALT_PNT']
zenith_mask = (zenith_min <= zenith) & (zenith < zenith_max)
print(repr(zenith_mask), zenith_mask.shape)
filtered_obs_table1 = obs_table1[zenith_mask]
print(repr(filtered_obs_table1))

# case 2: using `~gammapy.make.make_test_observation_table` (fails!)
#         this function returns a gammapy ObservationTable object
observatory_name='HESS'
n_obs = 10
obs_table2 = make_test_observation_table(observatory_name, n_obs)
zenith_min = Angle(20., 'degree')
zenith_max = Angle(40., 'degree')
zenith = Angle(90., 'degree') - obs_table2['ALT']
#zenith = Angle(90., 'degree') - Angle(obs_table2['ALT'], 'deg') # neede if returning obs_table as astropy Table (using array object: return Table(np.array(obs_table)))
zenith_mask = (zenith_min <= zenith) & (zenith < zenith_max)
print(repr(zenith_mask), zenith_mask.shape)
filtered_obs_table2 = obs_table2[zenith_mask]
print(repr(filtered_obs_table2))

import IPython; IPython.embed()

##EL PROBLEMA ESTA EN make_test_observation_table!!!!!!!!!!!!!
## en realidad no: con data_store.make_observation_table funciona porque devuelve una Table; pero make_test_observation_table devuelve una ObservationTable; si lo fuerzo a devolver una Table funciona
## asi que el problema esta en la clase ObservationTable!!!!!!
#TODO: hacer que ambas funciones (make_test_observation_table y data_store.make_observation_table) devuelvan una ObservationTable!!!
# hacer que ObservationTable funcione bien!!!!!!!!!!!!!!!!!!!!!

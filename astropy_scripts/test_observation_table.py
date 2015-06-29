from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
import numpy as np
from astropy.table import Table
from astropy.coordinates import Angle
from gammapy.utils.random import sample_sphere
from gammapy.obs import DataStore, ObservationTable
from gammapy.datasets.make import make_test_observation_table

# minimum (not?) working example :-)
# with observation table
print()
print("minimum example ObservationTable")
obs_table = ObservationTable()
obs_id = np.arange(10)
obs_table['OBS_ID'] = obs_id
obs_table.meta['NAME'] = 'table' # this line makes it fail!!! (really?) -> It works!!! :-)
mask = obs_table['OBS_ID'] == 1
print(repr(mask), mask.shape)
filtered_obs_table = obs_table[mask]
print(repr(filtered_obs_table))

# minimum working example :-)
# with astropy table
print()
print("minimum example Table")
table = Table()
obs_id = np.arange(10)
table['OBS_ID'] = obs_id
table.meta['NAME'] = 'table' # this line here works
mask = table['OBS_ID'] == 1
print(repr(mask), mask.shape)
filtered_table = table[mask]
print(repr(filtered_table))

# case 0: usin ObservationTable directly (works)
print()
print("case 0")
obs_table0 = ObservationTable()
n_obs_start = 1
n_obs = 10
obs_id = np.arange(n_obs_start, n_obs_start + n_obs)
obs_table0['OBS_ID'] = obs_id
az, alt = sample_sphere(len(obs_id), (0, 360), (45, 90), 'degree')
az = Angle(az, 'degree')
alt = Angle(alt, 'degree')
obs_table0['AZ'] = az
obs_table0['ALT'] = alt
zenith_min = Angle(20., 'degree')
zenith_max = Angle(40., 'degree')
zenith = Angle(90., 'degree') - obs_table0['ALT']
zenith_mask = (zenith_min <= zenith) & (zenith < zenith_max)
print(repr(zenith_mask), zenith_mask.shape)
filtered_obs_table0 = obs_table0[zenith_mask]
print(repr(filtered_obs_table0))

# case 1: using `~gammapy.obs.DataStore.make_observation_table` (works)
#         this function returns an astropy Table object
print()
print("case 1")
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

# case 2: using `~gammapy.make.make_test_observation_table` (fails!) -> fixed
#         this function returns a gammapy ObservationTable object
print()
print("case 2")
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

##import IPython; IPython.embed()

##EL PROBLEMA ESTA EN make_test_observation_table!!!!!!!!!!!!! -> fixed
## en realidad no: con data_store.make_observation_table funciona porque devuelve una Table; pero make_test_observation_table devuelve una ObservationTable; si lo fuerzo a devolver una Table funciona
## asi que el problema esta en la clase ObservationTable!!!!!!
#TODO: hacer que ambas funciones (make_test_observation_table y data_store.make_observation_table) devuelvan una ObservationTable!!!
# hacer que ObservationTable funcione bien!!!!!!!!!!!!!!!!!!!!!
# AL USAR DIRECTAMENTE UNA OBS_TABLE (simple) FUNIONA!!!! -> probar obs_table compleja; el problema podria estar en el make_test_observation_table !!! (intentar usar make_test_observation_table con una obs table sencilla?)!!! -> ya esta: el problema era un fallo tonto al definir el header de la tabla.

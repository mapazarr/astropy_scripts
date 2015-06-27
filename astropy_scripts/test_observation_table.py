#import gammapy
from gammapy.obs import DataStore
from astropy.coordinates import Angle
from gammapy.datasets.make import make_test_observation_table
HESSFITS_MPP = '/home/mapaz/astropy/gammapy_tutorial/HESS_fits_data/pa/Model_Deconvoluted_Prod26/Mpp_Std/'
data_store = DataStore(dir=HESSFITS_MPP)
observation_table = data_store.make_observation_table()
#observatory_name='HESS'
#n_obs = 100
#observation_table = make_test_observation_table(observatory_name, n_obs)
zenith_edges = Angle([0., 20., 40.], 'degree')
zenith = Angle(90., 'degree') - observation_table['ALT_PNT']
#zenith = Angle(90., 'degree') - observation_table['ALT']
for i in range (0, len(zenith_edges) - 1):
    zenith_mask = (zenith_edges[i] <= zenith) & (zenith < zenith_edges[i + 1])
    print(repr(zenith_mask), zenith_mask.shape)
    observation_table_filtered = observation_table[zenith_mask]
gammapy.obs

from gammapy.obs import ObservationTable
from gammapy.datasets.make import make_test_observation_table
from astropy.coordinates import Angle
observatory_name='HESS'
n_obs = 10
obs_table = make_test_observation_table(observatory_name, n_obs)
zenith_min = Angle(20., 'degree')
zenith_max = Angle(40., 'degree')
#zenith = Angle(90., 'degree') - obs_table['ALT']
zenith = Angle(90., 'degree') - Angle(obs_table['ALT'], 'deg')
zenith_mask = (zenith_min <= zenith) & (zenith < zenith_max)
print(repr(zenith_mask), zenith_mask.shape)
filtered_obs_table = obs_table[zenith_mask]

##EL PROBLEMA ESTA EN make_test_observation_table!!!!!!!!!!!!!
## en realidad no: con data_store.make_observation_table funciona porque devuelve una Table; pero make_test_observation_table devuelve una ObservationTable; si lo fuerzo a devolver una Table funciona
## asi que el problema esta en la clase ObservationTable!!!!!!
#TODO: hacer que ambas funciones (make_test_observation_table y data_store.make_observation_table) devuelvan una ObservationTable!!!
# hacer que ObservationTable funcione bien!!!!!!!!!!!!!!!!!!!!!

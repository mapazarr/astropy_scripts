from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
from astropy.table import Table
from gammapy.obs import (convert_obs_list_format_to_gammapy,
                         convert_obs_list_hess_to_gammapy)
from gammapy.obs import DataStore

# test: convert the H.E.S.S. runinfo.fits observation list
# from the PA FITS production

fits_path = '/home/mapaz/HESS/fits_data/pa_fits_prod02/pa/Model_Deconvoluted_Prod26/Mpp_Std/'
infile = fits_path + 'runinfo.fits'

obs_list = Table.read(infile)

print()
print("obs_list")
print(obs_list.meta)
print(obs_list)

obs_table = convert_obs_list_hess_to_gammapy(obs_list)

print()
print("obs_table")
print(obs_table.meta)
print(obs_table)

obs_table = convert_obs_list_format_to_gammapy(obs_list, scheme='hess')

print()
print("obs_table")
print(obs_table.meta)
print(obs_table)

# TODO: save a copy of the file? -> not necessary (since I am adapting the datastore class)?

# test how the datastore classes handle the converter

# if I try to import the converters in datastore.py, I get errors
# I think it is a circular dependency (I also can't import make_test_observation_table in there for instance
# moving the converters (convert.py) from datasets to obs didn't help
# I am updating the corresponding __init__.py accordingly...
# moving the functions to observation.py seems to work

data_store = DataStore(dir=fits_path)
obs_table = data_store.make_observation_table()
print()
print("obs_table datastore")
print(obs_table.meta)
print(obs_table)

event_list_files = data_store.make_table_of_files(obs_table, filetypes=['events'])
print()
print("event_list_files")
print(event_list_files)

aeff_list_files = data_store.make_table_of_files(obs_table, filetypes=['effective area'])
print()
print("aeff_list_files")
print(aeff_list_files)

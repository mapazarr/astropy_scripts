from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
from astropy.table import Table
from gammapy.datasets import (convert_obs_list_format_to_gammapy,
                              convert_obs_list_hess_to_gammapy)
##from gammapy.obs import DataStore

# test: convert the H.E.S.S. runinfo.fits observation list
# from the PA FITS production

inputdir = '/home/mapaz/HESS/fits_data/pa_fits_prod02/pa/Model_Deconvoluted_Prod26/Mpp_Std/'
infile = inputdir + 'runinfo.fits'

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

###data_store = DataStore(dir=fits_path)
###obs_table = data_store.make_observation_table()
###print()
###print("obs_table")
###print(obs_table.meta)
###print(obs_table)
###
###event_list_files = data_store.make_table_of_files(observation_table, filetypes=['events'])
###aeff_list_files = data_store.make_table_of_files(observation_table, filetypes=['effective area'])



#import IPython; IPython.embed()


# TODO: I probably have to "update" gammapy/obs/datastore.py to use the converter!!!
# DataStore.__init__ defines the file name
# DataStoreIndexTable.read is the one reading the file
# DataStoreIndexTable._init_cleanup -> review!!!
# check other TODOs within this PR (i.e. glon glat no units, etc)
#
# does runinfo have units??? or does the  DataStore/DataStoreIndexTable classes add them???!!!

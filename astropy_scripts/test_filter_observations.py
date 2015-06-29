from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
import numpy as np
from astropy.coordinates import Angle
import gammapy
from gammapy.obs import ObservationTable
from gammapy.datasets.make import make_test_observation_table

def test_filter_observations():
    # create random observation table
    observatory_name='HESS'
    n_obs = 10
    obs_table = make_test_observation_table(observatory_name, n_obs)
    print("obs_table")
    print(obs_table)

    # test apply a mask to the table (in zenith angle)
    print()
    print("Test apply mask in zenith angle:")
    zenith_min = Angle(20., 'degree')
    zenith_max = Angle(40., 'degree')
    zenith = Angle(90., 'degree') - obs_table['ALT']
    zenith_mask = (zenith_min <= zenith) & (zenith < zenith_max)
    filtered_obs_table = obs_table[zenith_mask]
    print("obs_table")
    print(obs_table)
    print("filtered_obs_table")
    print(filtered_obs_table)

    # test no selection: input and output tables should be the same
    print()
    print("Test no selection:")
    filtered_obs_table = obs_table.filter_observations()
    print("obs_table")
    print(obs_table)
    print("filtered_obs_table")
    print(filtered_obs_table)

#    return

    ##assert len(filtered_obs_table) == len(obs_table) #FIX THIS (IN CASE I DON'T RETURN A COPY OF THE OBJECT IN filter_observations!!!!!!!!!!!)!!!!!!!!!!
    #TODO: a lo mejor el fallo ya estaba de antes: intentar recuperar el estado inicial y hacer un test de select sky box sin mover el codigo!!!!!!!!!!!!!!!!!!!!!
    # intentar poner un ipithon en el table(mask)
    # intentar hacer table(mask=mask) en el sky box
    # probar en generico un table(mask) sin estar enterrado en 1000 clases
    # probar un obs_table(mask) tb.
    

    # filter some pars and check the correspoding values in the columns

    # test box selection in gal coordinates
    print()
    print("Test box selection in gal coordinates:")
    lon_min = -100.
    lon_max = 50.
    lat_min = -5.
    lat_max = 5.
    selection = dict(shape='box', frame='galactic',
                     lon=(lon_min, lon_max), lat=(lat_min, lat_max), border=2.)
    filtered_obs_table = obs_table.filter_observations(selection)
    print("obs_table")
    print(obs_table)
    print("filtered_obs_table")
    print(filtered_obs_table)


    # TODO: more!!! (implement some filters and test them!!!)

    # TODO: implement "inverse" filter: keep runs not inside, but outside selection criteria!!!

if __name__ == '__main__':
    test_filter_observations()

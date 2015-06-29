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
    assert len(filtered_obs_table) == len(obs_table)

    # filter some pars and check the correspoding values in the columns

    # test box selection in obs_id
    print()
    print("Test box selection in OBS_ID")
    variable = 'OBS_ID'
    min = 2
    max = 5
    selection = dict(shape='box', variable=variable,
                     min=min, max=max)
    filtered_obs_table = obs_table.filter_observations(selection)
    print("obs_table")
    print(obs_table)
    print("filtered_obs_table")
    print(filtered_obs_table)
    assert len(filtered_obs_table) == 3
    assert (min <= filtered_obs_table[variable]).all()
    assert (filtered_obs_table[variable] < max).all()

    # test box selection in obs_id inverted
    print()
    print("Test box selection in OBS_ID inverted")
    variable = 'OBS_ID'
    min = 2
    max = 5
    selection = dict(shape='box', variable=variable,
                     min=min, max=max, inverted=True)
    filtered_obs_table = obs_table.filter_observations(selection)
    print("obs_table")
    print(obs_table)
    print("filtered_obs_table")
    print(filtered_obs_table)
    assert len(filtered_obs_table) == 7
    assert ((min > filtered_obs_table[variable]) |
            (filtered_obs_table[variable] >= max)).all()

    # test circle selection in obs_id
    print()
    print("Test circle selection in OBS_ID")
    variable = 'OBS_ID'
    center = 4
    radius = 2
    selection = dict(shape='circle', variable=variable,
                     center=center, radius=radius)
    filtered_obs_table = obs_table.filter_observations(selection)
    print("obs_table")
    print(obs_table)
    print("filtered_obs_table")
    print(filtered_obs_table)
    assert len(filtered_obs_table) == 4
    assert (center - radius <= filtered_obs_table[variable]).all()
    assert (filtered_obs_table[variable] < center + radius).all()

    return

    # test sky box selection in gal coordinates
    print()
    print("Test sky box selection in gal coordinates:")
    lon_min = -100.
    lon_max = 50.
    lat_min = -5.
    lat_max = 5.
    selection = dict(shape='sky_box', frame='galactic',
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

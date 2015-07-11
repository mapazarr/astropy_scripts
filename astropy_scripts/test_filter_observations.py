from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
import numpy as np
from astropy.units import Quantity
from astropy.coordinates import Angle
from astropy.time import Time
from gammapy.obs import ObservationTable
from gammapy.datasets.make import make_test_observation_table
from gammapy.time import absolute_time

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
    zenith_max = Angle(30., 'degree')
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
    value_min = 2
    value_max = 5
    selection = dict(shape='box', variable=variable,
                     value_min=value_min, value_max=value_max)
    filtered_obs_table = obs_table.filter_observations(selection)
    print("obs_table")
    print(obs_table)
    print("filtered_obs_table")
    print(filtered_obs_table)
    assert len(filtered_obs_table) == 2
    assert (value_min < filtered_obs_table[variable]).all()
    assert (filtered_obs_table[variable] < value_max).all()

    # test box selection in obs_id inverted
    print()
    print("Test box selection in OBS_ID inverted")
    #variable = 'OBS_ID'
    #value_min = 2
    #value_max = 5
    selection = dict(shape='box', variable=variable,
                     value_min=value_min, value_max=value_max, inverted=True)
    filtered_obs_table = obs_table.filter_observations(selection)
    print("obs_table")
    print(obs_table)
    print("filtered_obs_table")
    print(filtered_obs_table)
    assert len(filtered_obs_table) == 8
    assert ((value_min >= filtered_obs_table[variable]) |
            (filtered_obs_table[variable] >= value_max)).all()

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
    assert len(filtered_obs_table) == 3
    assert (center - radius < filtered_obs_table[variable]).all()
    assert (filtered_obs_table[variable] < center + radius).all()

    # test box selection in alt
    print()
    print("Test box selection in ALT")
    variable = 'ALT'
    #value_min = 60.
    #value_max = 70.
    value_min = Angle(60., 'degree')
    value_max = Angle(70., 'degree')
    #value_min = Angle(60., 'degree').to('radian')
    #value_max = Angle(70., 'degree').to('radian')
    selection = dict(shape='box', variable=variable,
                     value_min=value_min, value_max=value_max)
    filtered_obs_table = obs_table.filter_observations(selection)
    print("obs_table")
    print(obs_table)
    print("filtered_obs_table")
    print(filtered_obs_table)
    #assert (value_min < Quantity(filtered_obs_table[variable])).all()
    #assert (Quantity(filtered_obs_table[variable]) < value_max).all()
    assert (value_min < Angle(filtered_obs_table[variable])).all()
    assert (Angle(filtered_obs_table[variable]) < value_max).all()

    # test box selection in zenith angle
    print()
    print("Test box selection in zenith")
    variable = 'zenith'
    value_min = Angle(20., 'degree')
    value_max = Angle(30., 'degree')
    selection = dict(shape='box', variable=variable,
                     value_min=value_min, value_max=value_max)
    filtered_obs_table = obs_table.filter_observations(selection)
    print("obs_table")
    print(obs_table)
    print("filtered_obs_table")
    print(filtered_obs_table)
    zenith = Angle(90., 'degree') - filtered_obs_table['ALT']
    print("zenit = {}".format(zenith))
    assert (value_min < zenith).all()
    assert (zenith < value_max).all()

    # test box selection in time_start
    print()
    print("Test box selection in TIME_START")
    variable = 'TIME_START'
    value_min = Time('2012-01-01 00:00:00', format='iso', scale='utc')
    value_max = Time('2014-01-01 00:00:00', format='iso', scale='utc')
    selection = dict(shape='box', variable=variable,
                     value_min=value_min, value_max=value_max)
    filtered_obs_table = obs_table.filter_observations(selection)
    print("obs_table")
    print(obs_table)
    print("filtered_obs_table")
    print(filtered_obs_table)
    time_start = absolute_time(filtered_obs_table['TIME_START'], filtered_obs_table.meta)
    print("time_start = {}".format(repr(time_start)))
    assert (value_min < time_start).all()
    assert (time_start < value_max).all()

    # test box selection in time: (time_start, time_stop) within (value_min, value_min)
    print()
    print("Test box selection in time")
    # new obs table with very close (in time) observations (and times in absolute times)
    datestart = Time('2012-01-01 00:03:00', format='iso', scale='utc')
    dateend = Time('2012-01-01 02:03:00', format='iso', scale='utc')
    obs_table_time = make_test_observation_table(observatory_name, n_obs,
                                            datestart, dateend, True)
    variable = 'time'
    value_min = Time('2012-01-01 01:00:00', format='iso', scale='utc')
    value_max = Time('2012-01-01 02:00:00', format='iso', scale='utc')
    selection = dict(shape='box', variable=variable,
                     value_min=value_min, value_max=value_max)
    filtered_obs_table = obs_table_time.filter_observations(selection)
    print("obs_table_time")
    print(obs_table_time)
    print("filtered_obs_table")
    print(filtered_obs_table)
    if filtered_obs_table.meta['TIME_FORMAT'] == 'absolute': # this is true, since I asked for it when creating the obs table
        time_start = filtered_obs_table['TIME_START']
        time_stop = filtered_obs_table['TIME_STOP']
    else:
        time_start = absolute_time(filtered_obs_table['TIME_START'], filtered_obs_table.meta)
        time_stop = absolute_time(filtered_obs_table['TIME_STOP'], filtered_obs_table.meta)
    print("time_start = {}".format(repr(time_start)))
    print("time_stop = {}".format(repr(time_stop)))
    assert (value_min < time_start).all()
    assert (time_start < value_max).all()
    assert (value_min < time_stop).all()
    assert (time_stop < value_max).all()

    return


    # TODO: continue with the sky_box/sky_circle tests!!!!


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





if __name__ == '__main__':
    test_filter_observations()

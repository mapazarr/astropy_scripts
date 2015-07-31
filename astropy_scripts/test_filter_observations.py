from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
import numpy as np
from astropy.units import Quantity
from astropy.coordinates import Angle, SkyCoord
from astropy.time import Time
from gammapy.obs import ObservationTable
from gammapy.datasets.make import make_test_observation_table
from gammapy.time import absolute_time
from gammapy.catalog import skycoord_from_table

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
    selection = dict(type='par_box', variable=variable,
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
    selection = dict(type='par_box', variable=variable,
                     value_min=value_min, value_max=value_max, inverted=True)
    filtered_obs_table = obs_table.filter_observations(selection)
    print("obs_table")
    print(obs_table)
    print("filtered_obs_table")
    print(filtered_obs_table)
    assert len(filtered_obs_table) == 8
    assert ((value_min >= filtered_obs_table[variable]) |
            (filtered_obs_table[variable] >= value_max)).all()

##    # test circle selection in obs_id
##    print()
##    print("Test circle selection in OBS_ID")
##    variable = 'OBS_ID'
##    center = 4
##    radius = 2
##    selection = dict(shape='circle', variable=variable,
##                     center=center, radius=radius)
##    filtered_obs_table = obs_table.filter_observations(selection)
##    print("obs_table")
##    print(obs_table)
##    print("filtered_obs_table")
##    print(filtered_obs_table)
##    assert len(filtered_obs_table) == 3
##    assert (center - radius < filtered_obs_table[variable]).all()
##    assert (filtered_obs_table[variable] < center + radius).all()

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
    selection = dict(type='par_box', variable=variable,
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

##    # test box selection in zenith angle
##    print()
##    print("Test box selection in zenith")
##    variable = 'zenith'
##    value_min = Angle(20., 'degree')
##    value_max = Angle(30., 'degree')
##    selection = dict(shape='box', variable=variable,
##                     value_min=value_min, value_max=value_max)
##    filtered_obs_table = obs_table.filter_observations(selection)
##    print("obs_table")
##    print(obs_table)
##    print("filtered_obs_table")
##    print(filtered_obs_table)
##    zenith = Angle(90., 'degree') - filtered_obs_table['ALT']
##    print("zenit = {}".format(zenith))
##    assert (value_min < zenith).all()
##    assert (zenith < value_max).all()

##    # test box selection in time_start
##    print()
##    print("Test box selection in TIME_START")
##    variable = 'TIME_START'
##    value_min = Time('2012-01-01 00:00:00', format='iso', scale='utc')
##    value_max = Time('2014-01-01 00:00:00', format='iso', scale='utc')
##    selection = dict(shape='box', variable=variable,
##                     value_min=value_min, value_max=value_max)
##    filtered_obs_table = obs_table.filter_observations(selection)
##    print("obs_table")
##    print(obs_table)
##    print("filtered_obs_table")
##    print(filtered_obs_table)
##    time_start = absolute_time(filtered_obs_table['TIME_START'], filtered_obs_table.meta)
##    print("time_start = {}".format(repr(time_start)))
##    assert (value_min < time_start).all()
##    assert (time_start < value_max).all()

    # test box selection in time: (time_start, time_stop) within (value_min, value_max)
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
    selection = dict(type='time_box', variable=variable,
                     time_min=value_min, time_max=value_max)
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

    # test sky box selection in gal coordinates
    print()
    print("Test sky box selection in gal coordinates:")
    lon_range = Angle([-100., 50.], 'degree')
    lat_range = Angle([-25., 25.], 'degree')
    frame = 'galactic'
    border = Angle(2., 'degree')
    selection = dict(type='sky_box', frame=frame,
                     lon=(lon_range[0], lon_range[1]),
                     lat=(lat_range[0], lat_range[1]),
                     border=border)
    lon_range_eff = (lon_range[0] - border, lon_range[1] + border)
    lat_range_eff = (lat_range[0] - border, lat_range[1] + border)
    print("Applied filters {0} {1}".format(lon_range_eff, lat_range_eff))
    filtered_obs_table = obs_table.filter_observations(selection)
    print("obs_table")
    print(obs_table)
    skycoord = skycoord_from_table(obs_table)
    skycoord = skycoord.transform_to(frame)
    print("Transformed coords: {}".format(skycoord))
    if any(l < 0 for l in lon_range_eff):
        print("Wrapping lon: {}".format(skycoord.data.lon.wrap_at(Angle(180, 'degree')).to('degree')))
    print("filtered_obs_table")
    print(filtered_obs_table)
    skycoord = skycoord_from_table(filtered_obs_table)
    skycoord = skycoord.transform_to(frame)
    print("Transformed filtered coords: {}".format(skycoord))
    if any(l < 0 for l in lon_range_eff):
        print("Wrapping lon: {}".format(skycoord.data.lon.wrap_at(Angle(180, 'degree')).to('degree')))
    lon = skycoord.data.lon
    lat = skycoord.data.lat
    if any(l < 0 for l in lon_range_eff):
        lon = lon.wrap_at(Angle(180, 'degree'))
    assert ((lon_range_eff[0] < lon) & (lon < lon_range_eff[1]) &
            (lat_range_eff[0] < lat) & (lat < lat_range_eff[1])).all()
    # test on the inverted selection
    selection['inverted'] = True
    print('selection', selection)
    inv_filtered_obs_table = obs_table.filter_observations(selection)
    print("inv_filtered_obs_table")
    print(inv_filtered_obs_table)
    skycoord = skycoord_from_table(inv_filtered_obs_table)
    skycoord = skycoord.transform_to(frame)
    lon = skycoord.data.lon
    lat = skycoord.data.lat
    if any(l < 0 for l in lon_range_eff):
        lon = lon.wrap_at(Angle(180, 'degree'))
    assert ((lon_range_eff[0] >= lon) | (lon >= lon_range_eff[1]) |
            (lat_range_eff[0] >= lat) | (lat >= lat_range_eff[1])).all()
    # the sum of number of entries in both selections should be the total number of entries
    print("sum = {}; original = {}".format(len(filtered_obs_table) + len(inv_filtered_obs_table), len(obs_table)))
    assert len(filtered_obs_table) + len(inv_filtered_obs_table) == len(obs_table)

    return




    # test sky box selection in radec coordinates
    print()
    print("Test sky box selection in radec coordinates:")
    lon_range = Angle([150., 300.], 'degree')
    lat_range = Angle([-50., 0.], 'degree')
    frame = 'icrs'
    border = Angle(2., 'degree')
    selection = dict(type='sky_box', frame=frame,
                     lon=(lon_range[0], lon_range[1]),
                     lat=(lat_range[0], lat_range[1]),
                     border=border)
    lon_range_eff = (lon_range[0] - border, lon_range[1] + border)
    lat_range_eff = (lat_range[0] - border, lat_range[1] + border)
    print("Applied filters {0} {1}".format(lon_range_eff, lat_range_eff))
    filtered_obs_table = obs_table.filter_observations(selection)
    print("obs_table")
    print(obs_table)
    skycoord = skycoord_from_table(obs_table)
    skycoord = skycoord.transform_to(frame)
    print("Transformed coords: {}".format(skycoord))
    print("filtered_obs_table")
    print(filtered_obs_table)
    skycoord = skycoord_from_table(filtered_obs_table)
    skycoord = skycoord.transform_to(frame)
    print("Transformed filtered coords: {}".format(skycoord))
    lon = skycoord.data.lon
    lat = skycoord.data.lat
    assert ((lon_range_eff[0] < lon) & (lon < lon_range_eff[1])).all()
    assert ((lat_range_eff[0] < lat) & (lat < lat_range_eff[1])).all()

    #####import IPython; IPython.embed()

    # TODO: test sky box selection of a box with lon on the range (0, 180) deg, with border 2 !!!!
    # also (0, 10) deg, with border 2 !!!!
    # also (0, 200) deg, with border 2 !!!!
    # what is border for???!!!

    # TODO: continue with the sky_circle tests!!!!
    # TODO: seed all tests with random obs tables (or random bg cubes?) otherwise I get strange failures from time to time!!!!!!

    # TODO: I'm testing on the elements that survive the selection, but not on the ones that fail!!!!
    #       if for instance all elements are skipped because of a problem in the selection, then the test passes as good!!!
    #       would testing on the inverted selection help??!!!
    #       not completely: I can have the same problem!!!
    #       test that the sum of entries in both tables should be the total number of entries in the original table

    # TODO: a function to test coordinate filters, then call it everywhere




def test_find_observations():
    '''
    gammapy-find-obs -h
    gammapy-find-obs runinfo.fits
    gammapy-find-obs runinfo.fits --x 3
    gammapy-find-obs runinfo.fits --x 3 --r 2
    gammapy-find-obs runinfo.fits --x 3 --y 4 --r 2 --system 'radec'
    gammapy-find-obs runinfo.fits --x 3 --y 4 --dx 5 --dy 6 --system 'radec'
    gammapy-find-obs runinfo.fits --t_start 2 --t_stop 3
    gammapy-find-obs runinfo.fits --t_start '2012-01-01 00:00:00' --t_stop '2014-01-01 00:00:00'
    gammapy-find-obs runinfo.fits --t_start '2012-01-01 01:00:00' --t_stop '2012-01-01 02:00:00' # needs "special" test obs table
    gammapy-find-obs runinfo.fits --par_name "OBS_ID"
    gammapy-find-obs runinfo.fits --par_name "OBS_ID" --par_min 2 --par_max 6
    gammapy-find-obs runinfo.fits --par_name "ALT" --par_min 60 --par_max 70
    '''

    # TODO: test the time selection using an obs table in 2 different formats (absolute and MET)!!!!! -> actually not necessary, it's tested in test_filter_observations
    # TODO: a really long selection, with many criteria


def test_skycoord_angle_wrapping():
    """Test astropy SkyCoord angle wrapping."""

    lon = Angle(0. - 50., 'degree')
    lat = Angle(0., 'degree')
    frame = 'galactic'
    sky_coord = SkyCoord(lon, lat, frame=frame)
    print(sky_coord)

    lon = Angle(0. - 50., 'degree')
    lat = Angle(0., 'degree')
    frame = 'icrs'
    sky_coord = SkyCoord(lon, lat, frame=frame)
    print(sky_coord)

    lon = Angle(360. - 50., 'degree')
    lat = Angle(0., 'degree')
    frame = 'galactic'
    sky_coord = SkyCoord(lon, lat, frame=frame)
    print(sky_coord)

    lon = Angle(360. - 50., 'degree')
    lat = Angle(0., 'degree')
    frame = 'icrs'
    sky_coord = SkyCoord(lon, lat, frame=frame)
    print(sky_coord)

    # It seems like SkyCoord always wraps at 360 deg.

    lon = sky_coord.data.lon
    lat = sky_coord.data.lat
    print(lon, lat)
    lon.wrap_at(Angle(180, 'degree'))
    print(lon, lat)
    lon = lon.wrap_at(Angle(180, 'degree'))
    print(lon, lat)


if __name__ == '__main__':
    test_filter_observations()
    test_find_observations()
    test_skycoord_angle_wrapping()

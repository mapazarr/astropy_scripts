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


def common_sky_region_select_test_routines(obs_table, selection):
    """Common routines for the tests of sky_box/sky_circle selection of obs tables."""
    type = selection['type']
    if type not in ['sky_box', 'sky_circle']:
        raise ValueError("Invalid type: {}".format(type))

    if type == 'sky_box':
        lon_range_eff = (selection['lon'][0] - selection['border'],
                         selection['lon'][1] + selection['border'])
        lat_range_eff = (selection['lat'][0] - selection['border'],
                         selection['lat'][1] + selection['border'])
        print("Applied selections {0} {1}".format(lon_range_eff, lat_range_eff))
    elif type == 'sky_circle':
        lon_cen = selection['lon']
        lat_cen = selection['lat']
        center = SkyCoord(lon_cen, lat_cen, frame=selection['frame'])
        radius_eff = selection['radius'] + selection['border']
        print("Applied selections {0} around {1})".format(radius_eff, center))

    do_wrapping = False
    # not needed in the case of sky_circle
    if (type == 'sky_box' and
        any(l < Angle(0., 'degree') for l in lon_range_eff)):
        do_wrapping = True

    # observation table
    print()
    print("obs_table")
    print('header:', obs_table.meta)
    print(obs_table)
    skycoord = skycoord_from_table(obs_table)
    if type == 'sky_box':
        skycoord = skycoord.transform_to(selection['frame'])
        print("Transformed coords: {}".format(skycoord))
        if do_wrapping:
            print("Wrapping lon: {}".format(skycoord.data.lon.wrap_at(Angle(180, 'degree')).to('degree')))
    elif type == 'sky_circle':
        ang_distance = skycoord.separation(center)
        print("Angular distances: {}".format(ang_distance))

    # test on the selection
    print()
    print("selected_obs_table")
    print('selection:', selection)
    selected_obs_table = obs_table.select_observations(selection)
    print(selected_obs_table)
    skycoord = skycoord_from_table(selected_obs_table)
    if type == 'sky_box':
        skycoord = skycoord.transform_to(selection['frame'])
        print("Transformed selected coords: {}".format(skycoord))
        if do_wrapping:
            print("Wrapping lon: {}".format(skycoord.data.lon.wrap_at(Angle(180, 'degree')).to('degree')))
        lon = skycoord.data.lon
        lat = skycoord.data.lat
        if do_wrapping:
            lon = lon.wrap_at(Angle(180, 'degree'))
        assert ((lon_range_eff[0] < lon) & (lon < lon_range_eff[1]) &
                (lat_range_eff[0] < lat) & (lat < lat_range_eff[1])).all()
    elif type == 'sky_circle':
        ang_distance = skycoord.separation(center)
        print("Angular selected distances: {}".format(ang_distance))
        assert (ang_distance < radius_eff).all()

    # test on the inverted selection
    print()
    print("inv_selected_obs_table")
    selection['inverted'] = True
    print('selection:', selection)
    inv_selected_obs_table = obs_table.select_observations(selection)
    print(inv_selected_obs_table)
    skycoord = skycoord_from_table(inv_selected_obs_table)
    if type == 'sky_box':
        skycoord = skycoord.transform_to(selection['frame'])
        print("Transformed inv selected coords: {}".format(skycoord))
        if do_wrapping:
            print("Wrapping lon: {}".format(skycoord.data.lon.wrap_at(Angle(180, 'degree')).to('degree')))
        lon = skycoord.data.lon
        lat = skycoord.data.lat
        if do_wrapping:
            lon = lon.wrap_at(Angle(180, 'degree'))
        assert ((lon_range_eff[0] >= lon) | (lon >= lon_range_eff[1]) |
                (lat_range_eff[0] >= lat) | (lat >= lat_range_eff[1])).all()
    elif type == 'sky_circle':
        ang_distance = skycoord.separation(center)
        print("Angular inv selected distances: {}".format(ang_distance))
        assert (ang_distance >= radius_eff).all()

    # the sum of number of entries in both selections should be the total number of entries
    print()
    print(" original = {}; selected = {}; inv_selected = {}; sum = {}"
          .format(len(obs_table), len(selected_obs_table), len(inv_selected_obs_table),
                  len(selected_obs_table) + len(inv_selected_obs_table)))
    assert len(selected_obs_table) + len(inv_selected_obs_table) == len(obs_table)


def test_select_observations():
    test_select_parameter_box()
    test_select_time_box()
    test_select_sky_regions()

def test_select_parameter_box():
    """Test select parameter box."""

    # create random observation table
    obs_table = make_test_observation_table(n_obs=10)
    print("obs_table")
    print(obs_table)

    # test apply a mask to the table (in zenith angle)
    print()
    print("Test apply mask in zenith angle:")
    zenith_min = Angle(20., 'degree')
    zenith_max = Angle(30., 'degree')
    zenith = Angle(90., 'degree') - obs_table['ALT']
    zenith_mask = (zenith_min <= zenith) & (zenith < zenith_max)
    selected_obs_table = obs_table[zenith_mask]
    print("obs_table")
    print(obs_table)
    print("selected_obs_table")
    print(selected_obs_table)

    # test no selection: input and output tables should be the same
    print()
    print("Test no selection:")
    selected_obs_table = obs_table.select_observations()
    print("obs_table")
    print(obs_table)
    print("selected_obs_table")
    print(selected_obs_table)
    assert len(selected_obs_table) == len(obs_table)

    # select some pars and check the correspoding values in the columns

    # test box selection in obs_id
    print()
    print("Test box selection in OBS_ID")
    variable = 'OBS_ID'
    value_range = [2, 5]
    selection = dict(type='par_box', variable=variable, value_range=value_range)
    selected_obs_table = obs_table.select_observations(selection)
    print("obs_table")
    print(obs_table)
    print("selected_obs_table")
    print(selected_obs_table)
    assert len(selected_obs_table) == 3
    assert (value_range[0] <= selected_obs_table[variable]).all()
    assert (selected_obs_table[variable] < value_range[1]).all()

    # test box selection in obs_id inverted
    print()
    print("Test box selection in OBS_ID inverted")
    #variable = 'OBS_ID'
    #value_min = 2
    #value_max = 5
    selection = dict(type='par_box', variable=variable,
                     value_range=value_range, inverted=True)
    selected_obs_table = obs_table.select_observations(selection)
    print("obs_table")
    print(obs_table)
    print("selected_obs_table")
    print(selected_obs_table)
    assert len(selected_obs_table) == 7
    assert ((value_range[0] > selected_obs_table[variable]) |
            (selected_obs_table[variable] >= value_range[1])).all()

    # test box selection in alt
    print()
    print("Test box selection in ALT")
    variable = 'ALT'
    #value_range = [60., 70.]
    value_range = Angle([60., 70.], 'degree')
    #value_range = Angle([60., 70.], 'degree').to('radian')
    selection = dict(type='par_box', variable=variable, value_range=value_range)
    selected_obs_table = obs_table.select_observations(selection)
    print("obs_table")
    print(obs_table)
    print("selected_obs_table")
    print(selected_obs_table)
    #assert (value_range[0] < Quantity(selected_obs_table[variable])).all()
    #assert (Quantity(selected_obs_table[variable]) < value_range[1]).all()
    assert (value_range[0] < Angle(selected_obs_table[variable])).all()
    assert (Angle(selected_obs_table[variable]) < value_range[1]).all()


def test_select_time_box():
    """Test select time box."""

    # create random observation table with very close (in time)
    # observations (and times in absolute times)
    datestart = Time('2012-01-01T00:30:00', format='isot', scale='utc')
    dateend = Time('2012-01-01T02:30:00', format='isot', scale='utc')
    obs_table_time = make_test_observation_table(n_obs=10,
                                                 datestart=datestart,
                                                 dateend=dateend,
                                                 use_abs_time=True)

    # test box selection in time: (time_start, time_stop) within value_range
    print()
    print("Test box selection in time")
    value_range = Time(['2012-01-01T01:00:00', '2012-01-01T02:00:00'],
                       format='isot', scale='utc')
    selection = dict(type='time_box', time_range=value_range)
    selected_obs_table = obs_table_time.select_observations(selection)
    print("obs_table_time")
    print(obs_table_time)
    print("selected_obs_table")
    print(selected_obs_table)
    if selected_obs_table.meta['TIME_FORMAT'] == 'absolute': # this is true, since I asked for it when creating the obs table
        time_start = selected_obs_table['TIME_START']
        time_stop = selected_obs_table['TIME_STOP']
    else:
        time_start = absolute_time(selected_obs_table['TIME_START'], selected_obs_table.meta)
        time_stop = absolute_time(selected_obs_table['TIME_STOP'], selected_obs_table.meta)
    print("time_start = {}".format(repr(time_start)))
    print("time_stop = {}".format(repr(time_stop)))
    assert (value_range[0] < time_start).all()
    assert (time_start < value_range[1]).all()
    assert (value_range[0] < time_stop).all()
    assert (time_stop < value_range[1]).all()


def test_select_sky_regions():
    """Test select sky regions (box and circle)."""

    # create random observation table with many entries
    obs_table = make_test_observation_table(n_obs=100)
    print("obs_table")
    print(obs_table)

    # test sky box selection in gal coordinates
    print()
    print("Test sky box selection in gal coordinates:")
    lon_range = Angle([-100., 50.], 'degree')
    lat_range = Angle([-25., 25.], 'degree')
    frame = 'galactic'
    border = Angle(2., 'degree')
    selection = dict(type='sky_box', frame=frame,
                     lon=lon_range,
                     lat=lat_range,
                     border=border)
    common_sky_region_select_test_routines(obs_table, selection)

    # test sky box selection in radec coordinates
    print()
    print("Test sky box selection in radec coordinates:")
    lon_range = Angle([150., 300.], 'degree')
    lat_range = Angle([-50., 0.], 'degree')
    frame = 'icrs'
    border = Angle(2., 'degree')
    selection = dict(type='sky_box', frame=frame,
                     lon=lon_range,
                     lat=lat_range,
                     border=border)
    common_sky_region_select_test_routines(obs_table, selection)

    # test sky circle selection in gal coordinates
    print()
    print("Test sky circle selection in gal coordinates:")
    lon_cen = Angle(0., 'degree')
    lat_cen = Angle(0., 'degree')
    radius = Angle(50., 'degree')
    frame = 'galactic'
    border = Angle(2., 'degree')
    selection = dict(type='sky_circle', frame=frame,
                     lon=lon_cen, lat=lat_cen,
                     radius=radius, border=border)
    common_sky_region_select_test_routines(obs_table, selection)

    # test sky circle selection in radec coordinates
    print()
    print("Test sky circle selection in gal coordinates:")
    lon_cen = Angle(130., 'degree')
    lat_cen = Angle(-40., 'degree')
    radius = Angle(50., 'degree')
    frame = 'icrs'
    border = Angle(2., 'degree')
    selection = dict(type='sky_circle', frame=frame,
                     lon=lon_cen, lat=lat_cen,
                     radius=radius, border=border)
    common_sky_region_select_test_routines(obs_table, selection)

    # TODO: test sky box selection of a box with lon on the range (0, 180) deg, with border 2 !!!!
    # also (0, 10) deg, with border 2 !!!!
    # also (0, 200) deg, with border 2 !!!!
    # what is border for???!!!

    # TODO: I'm testing on the elements that survive the selection, but not on the ones that fail!!!! -> ~done?
    #       if for instance all elements are skipped because of a problem in the selection, then the test passes as good!!!
    #       would testing on the inverted selection help??!!!
    #       not completely: I can have the same problem!!!
    #       on top of that, test that the sum of entries in both tables should be the total number of entries in the original table -> works, but is it enough?


def test_find_observations():
    """
    gammapy-find-obs -h
    gammapy-find-obs runinfo.fits
    gammapy-find-obs dummy.fits out.dummy.fits
    #gammapy-find-obs runinfo.fits --x 3
    #gammapy-find-obs runinfo.fits --x 3 --r 2
    #gammapy-find-obs runinfo.fits --x 3 --y 4 --r 2 --system 'icrs'
    gammapy-find-obs runinfo.fits --x 130 --y -40 --r 50 --system 'icrs'
    gammapy-find-obs runinfo.fits --x 0 --y 0 --r 50 --system 'galactic'
    #gammapy-find-obs runinfo.fits --x 3 --y 4 --dx 5 --dy 6 --system 'icrs'
    gammapy-find-obs runinfo.fits --x 225 --y -25 --dx 75 --dy 25 --system 'icrs'
    gammapy-find-obs runinfo.fits --x -25 --y 0 --dx 75 --dy 25 --system 'galactic'
    #gammapy-find-obs runinfo.fits --t_start 2 --t_stop 3
    gammapy-find-obs runinfo.fits --t_start '2012-01-01T00:00:00' --t_stop '2014-01-01T00:00:00'
    #gammapy-find-obs runinfo.fits --t_start '2012-01-01T01:00:00' --t_stop '2012-01-01T02:00:00' # needs "special" test obs table
    #gammapy-find-obs runinfo.fits --par_name 'OBS_ID'
    gammapy-find-obs runinfo.fits --par_name 'OBS_ID' --par_min 2 --par_max 6
    gammapy-find-obs runinfo.fits --par_name 'ALT' --par_min 60 --par_max 70
    """

    # TODO: a really long selection, with many criteria


def test_skycoord_angle_wrapping():
    """Test astropy SkyCoord angle wrapping."""

    # test SkyCoord wrapping

    print()

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

    # test Angle, Longitude, Latitude wrapping

    print()
    lon = sky_coord.data.lon
    lat = sky_coord.data.lat
    print(repr(lon), repr(lat))
    lon.wrap_at(Angle(180, 'degree')) # returns Angle
    print(repr(lon), repr(lat))
    lon = lon.wrap_at(Angle(180, 'degree')) # returns Angle
    print(repr(lon), repr(lat))

    print()
    lon = sky_coord.data.lon
    lat = sky_coord.data.lat
    print(repr(lon), repr(lat))
    lon.wrap_angle = Angle(180, 'degree') # updates Longitude
    print(repr(lon), repr(lat))

    # test angular distances

    print()
    sky_coord1 = SkyCoord(0., 0., unit='degree', frame='icrs')
    sky_coord2 = SkyCoord(0., 0., unit='degree', frame='icrs')
    print(sky_coord1.separation(sky_coord2))
    sky_coord1 = SkyCoord(5., 0., unit='degree', frame='icrs')
    sky_coord2 = SkyCoord(355., 0., unit='degree', frame='icrs')
    print(sky_coord1.separation(sky_coord2))
    sky_coord1 = SkyCoord(175., 0., unit='degree', frame='icrs')
    sky_coord2 = SkyCoord(-175., 0., unit='degree', frame='icrs')
    print(sky_coord1.separation(sky_coord2))
    sky_coord1 = SkyCoord(355., 0., unit='degree', frame='icrs')
    sky_coord2 = SkyCoord(-175., 0., unit='degree', frame='icrs')
    print(sky_coord1.separation(sky_coord2))
    sky_coord1 = SkyCoord(355., 0., unit='degree', frame='icrs')
    sky_coord2 = SkyCoord(-5., 0., unit='degree', frame='icrs')
    print(sky_coord1.separation(sky_coord2))
    sky_coord1 = SkyCoord(0., 85., unit='degree', frame='icrs')
    sky_coord2 = SkyCoord(0., -85., unit='degree', frame='icrs')
    print(sky_coord1.separation(sky_coord2))

    # distances work fine, no need to wrap


if __name__ == '__main__':
    test_select_observations()
    test_find_observations()
    test_skycoord_angle_wrapping()

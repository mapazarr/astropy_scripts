from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
import numpy as np
from astropy.coordinates import Angle, Longitude, Latitude
from astropy.tests.helper import assert_quantity_allclose
from gammapy.utils.random import sample_sphere


# test general case
np.random.seed(0)
lon, lat = sample_sphere(size=2)
assert_quantity_allclose(lon, Longitude([3.44829694, 4.49366732], 'radian'))
assert_quantity_allclose(lat, Latitude([0.20700192, 0.08988736], 'radian'))

#define test of setting lon_range 0, 360 and -180, 180

# test lon range explicitly (0, 360) deg
epsilon = 1.e-8
lon_range = Longitude([0., 360.-epsilon], 'degree')
#lon_range = Longitude([0., 360.], 'degree')
lon, lat = sample_sphere(size=100, lon_range=lon_range)
angle_0 = Angle(0., 'degree')
angle_360 = Angle(360., 'degree')
angle_m90 = Angle(-90., 'degree')
angle_90 = Angle(90., 'degree')
# test values in the desired range
assert ((angle_0 <= lon) & (lon < angle_360)).all()
assert ((angle_m90 <= lat) & (lat < angle_90)).all()
# test if values are distributed along the whole range
nbins = 4
angle_delta_0_360 = (angle_360 - angle_0)/nbins
angle_delta_m90_90 = (angle_90 - angle_m90)/nbins
for i in np.arange(nbins):
    assert ((angle_0 + i*angle_delta_0_360 <= lon) &
            (lon < angle_0 + (i + 1)*angle_delta_0_360)).any()
    assert ((angle_m90 + i*angle_delta_m90_90 <= lat) &
            (lat < angle_m90 + (i + 1)*angle_delta_m90_90)).any()

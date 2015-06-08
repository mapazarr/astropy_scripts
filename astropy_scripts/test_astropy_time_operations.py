from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
import numpy as np
from astropy.time import Time, TimeDelta
from gammapy.utils.time import time_ref_from_dict, time_relative_to_ref
from print_variable_debug import print_variable_debug as var_debug

# inspired from gammapy.utils.time

# define reference time in a dictionary
mjd_int = 500
mjd_frac = 0.5

time_ref_dict = dict(MJDREFI=mjd_int, MJDREFF=mjd_frac)

print("time_ref_dict")
print(time_ref_dict)
print(type(time_ref_dict))

print()

time_ref = time_ref_from_dict(time_ref_dict)

print("time_ref")
print(time_ref)
print(type(time_ref))

# test: is the time ref equal to the sum of mjd_int + mjd_frac?
decimal = 4
s_error = "time reference not compatible with defined values"
np.testing.assert_almost_equal(time_ref.mjd, mjd_int + mjd_frac, decimal, s_error)

print()

# define a time as time_ref + 1 sec
delta_time_1sec = TimeDelta(1., format='sec')

print("delta_time_1sec")
print(delta_time_1sec)
print(type(delta_time_1sec))

print()

time = time_ref + delta_time_1sec

print("time")
print(time)
print(type(time))

print()

delta_time_calc_per_hand = time - time_ref

print("delta_time_calc_per_hand")
print(delta_time_calc_per_hand)
print(type(delta_time_calc_per_hand))

print()

delta_time_calc_per_hand_sec = delta_time_calc_per_hand.sec

print("delta_time_calc_per_hand_sec")
print(delta_time_calc_per_hand_sec)
print(type(delta_time_calc_per_hand_sec))

print()

delta_time_calc_per_hand_format_sec = TimeDelta(time - time_ref, format='sec')

print("delta_time_calc_per_hand_format_sec")
print(delta_time_calc_per_hand_format_sec)
print(type(delta_time_calc_per_hand_format_sec))

# test if delta time is calculated correctly
delta_time = time_relative_to_ref(time, time_ref_dict)
decimal = 4
s_error = "delta time not compatible with defined value"
np.testing.assert_almost_equal(delta_time.sec, delta_time_1sec.sec, decimal, s_error)


# inspired from gammapy.datasets.make.make_test_observation_table()

print()

#build a reference
dateref = Time('2010-01-01 00:00:00', format='iso', scale='utc')

print("dateref")
var_debug(dateref)

print()

dateref_mjd = dateref.mjd

print("dateref_mjd")
var_debug(dateref_mjd)

print()

dateref_mjd_fra, dateref_mjd_int = np.modf(dateref_mjd)
print("dateref_mjd_int =", dateref_mjd_int)
print("dateref_mjd_fra =", dateref_mjd_fra)

print()

#header as dictionary
header = {'MJDREFI': dateref_mjd_int, 'MJDREFF': dateref_mjd_fra}

print("header")
var_debug(header)

print()

#convert a random date to seconds after the reference
date = Time('2011-11-11 11:11:11', format='iso', scale='utc')

print("date")
var_debug(date)

print()

date_rel_to_ref = time_relative_to_ref(date, header)

print("date_rel_to_ref")
var_debug(date_rel_to_ref)

print()

#build a delta time for the night of 6 h (from 22 h to 4 h)
delta_time_night_duration = TimeDelta(6.*60.*60., format='sec')

print("delta_time_night_duration")
var_debug(delta_time_night_duration)

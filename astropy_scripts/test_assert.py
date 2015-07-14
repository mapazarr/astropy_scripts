# Testing assert_allclose (copy/paste code in ipython)
#
# In case of assertion errors on arrays of quantity objects, such
# as `astropy.units.Quantity` or `~astropy.coordinates.Angle`,
# `numpy.testing.assert_allclose` will not show the values. It
# shows `[repr failed]`.
#
# For debugging you can use the `.value` methods of these object
# types. Don't forget to remove the call to `.value` for final
# code.

from numpy.testing import assert_allclose
from astropy.coordinates import Angle
from astropy.units import Quantity

a = Angle(0., 'degree')
b = Angle(1., 'degree')
assert_allclose(a, b) # ok
a2 = Angle([0., 0.], 'degree')
b2 = Angle([1., 1.], 'degree')
assert_allclose(a2, b2) # problem
q = Quantity(0., 'm')
r = Quantity(1., 'm')
assert_allclose(q, r) # ok
q2 = Quantity([0., 0.], 'm')
r2 = Quantity([1., 1.], 'm')
assert_allclose(q2, r2) # problem
#history

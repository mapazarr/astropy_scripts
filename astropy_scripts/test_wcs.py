from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
from astropy import wcs
from astropy.wcs import WCS
from astropy.coordinates import Angle
from astropy.tests.helper import assert_quantity_allclose
from gammapy.utils.wcs import (linear_wcs_to_arrays,
                               linear_arrays_to_wcs)

# ASTROPY STUFF

# test defining wcs object as astropy.wcs

w_wcs = wcs.WCS(naxis=2)
w_wcs.wcs.ctype = ["X", "Y"]
w_wcs.wcs.cunit = ["deg", "deg"]
w_wcs.wcs.crpix = [0.5, 0.5]
w_wcs.wcs.crval = [1.23, 4.56]
w_wcs.wcs.cdelt = [1.5, 1.5]

print("w_wcs")
print(w_wcs)
w_wcs.printwcs() # no units printed?
w_wcs_header = w_wcs.to_header()
print("w_wcs_header")
print(repr(w_wcs_header))

print()

# test defining wcs object as astropy.wcs.WCS

w_WCS = WCS(naxis=2)
# keywords appears then empty
#w_WCS.ctype = ["X", "Y"]
#w_WCS.cunit = ["deg", "deg"]
#w_WCS.crpix = [0.5, 0.5]
#w_WCS.crval = [1.23, 4.56]
#w_WCS.cdelt = [1.5, 1.5]
# keywords are filled
w_WCS.wcs.ctype = ["X", "Y"]
w_WCS.wcs.cunit = ["deg", "deg"]
w_WCS.wcs.crpix = [0.5, 0.5]
w_WCS.wcs.crval = [1.23, 4.56]
w_WCS.wcs.cdelt = [1.5, 1.5]

print("w_WCS")
print(w_WCS)
w_WCS.printwcs() # no units printed?
w_WCS_header = w_WCS.to_header()
print("w_WCS_header")
print(repr(w_WCS_header))

# GAMMAPY STUFF

# test WCS object
print()
print("Testing make WCS")
bins_x = Angle([1., 2., 3., 4.], 'degree')
#bins_y = Angle([-1.5, 0., 1.5], 'radian')
bins_y = Angle([-1.5, 0., 1.5], 'degree')
wcs = linear_arrays_to_wcs("X", "Y", bins_x, bins_y)
wcs_header = wcs.to_header()
print(repr(wcs_header))
print()
print("Testing recover bins")
nbins_x = len(bins_x) - 1
nbins_y = len(bins_y) - 1
reco_bins_x, reco_bins_y = linear_wcs_to_arrays(wcs, nbins_x, nbins_y)
print(repr(reco_bins_x))
print(repr(reco_bins_y))
# to deeply test the translation of the bin reconstruction, use this in the wcs creation function:
#w.wcs.crpix = [10, 15] # smthg different from origin (0.5, 0.5) !!!
# test: reconstructed bins should match original bins
assert_quantity_allclose(reco_bins_x, bins_x)
assert_quantity_allclose(reco_bins_y, bins_y)

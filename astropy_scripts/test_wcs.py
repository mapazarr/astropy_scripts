from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3

from astropy import wcs
from astropy.wcs import WCS

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

#import IPython; IPython.embed()

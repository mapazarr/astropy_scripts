
from gammapy.background import fill_acceptance_image
from gammapy.image import make_empty_image

from astropy.io import fits
from astropy.coordinates import SkyCoord, Angle
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
import astropy.units as u
from astropy.units.quantity import Quantity

import numpy as np
from scipy.integrate import quad
from matplotlib import pyplot as plt

def radial_gaussian_2D(r, sig):
    """Radially symmetric Gaussian function for emulating acceptance curve.

    Please give parameters in radians.

    TODO: I can't define parameters as `~astropy.coordinates.Angle` because it conflicts with the integral call afterwards!!!

    Parameters
    ----------
    r : `~numpy.ndarray`
    Radial coordinate.
    sig : float
    Gaussian width.

    Returns
    -------
    val: `~numpy.ndarray`
    Gaussian evaluated at the requested radii.
    """
    norm = 1.
    val = norm * np.exp(-np.power(r, 2.) / (2. * np.power(sig, 2.)))
    return val

#create empty image
image = make_empty_image()

#image is a astropy.io.fits.hdu.image.ImageHDU
# An ImageHDU has two important attributes:
#  - data, which behaves like a Numpy array, can be used to access the data
#  - header, which behaves like a dictionary, can be used to access the header information

print ""

#header returns a dictionary object
print "image.header"
#print image.header ##exactly as in fits file; it's equivalent to print image.header.tostring()
print image.header.tostring('\n', False, False)

print ""

#data returns a numpy array object
#print "image.data %s" %np.array_str(image.data)
print "image.data"
print(image.data)

print ""

#shape returns a tuple object
print "image.data.shape:", image.data.shape

print ""

#plot (empty) image
print "here comes the plot"
fig = plt.figure()
ax = fig.add_subplot(111)

ax.imshow(image.data[:,:], origin='lower')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('empty image')

plt.draw() #draw plot

print ""

#modify header entry
image.header['telescop'] = 'Python Observatory'
print image.header.tostring('\n', False, False)

print ""

#define center coordinate of the image
#in FITS the reference pixel is given by CRPIXn, and its coordinate by CRVALn
#the coordinate increment along the axis is given by CDELTn
x = image.header['CRVAL1']
y = image.header['CRVAL2']

#sanity checks: am I using galactic coordinates in degrees?
assert image.header['CTYPE1'] == 'GLON-CAR', "Error: x coordinate is not in Galactic coordinates!"
assert image.header['CTYPE2'] == 'GLAT-CAR', "Error: y coordinate is not in Galactic coordinates!"
assert image.header['CUNIT1'] == 'deg', "Error: x coordinate is not in degrees!"
assert image.header['CUNIT2'] == 'deg', "Error: y coordinate is not in degrees!"

center = SkyCoord(l=x*u.degree, b=y*u.degree, frame='galactic')

print "center"
print(center)

#define pixel sizes
x_pix_size = Angle(abs(image.header['CDELT1'])*u.degree)
y_pix_size = Angle(abs(image.header['CDELT2'])*u.degree)

print ""

#define radial acceptance and offset angles
offset = Angle(np.arange(0., 30., 0.1), unit=u.degree)
acceptance = np.zeros_like(offset)
sigma = Angle(1.0, unit=u.degree) #gaussian width
acceptance = radial_gaussian_2D(offset.to(u.radian).value, sigma.to(u.radian).value)

print "offset"
print(offset)
print "acceptance"
print(acceptance)

#fill acceptance in the image
#alternative: create a new object, instead of overwriting the empty image
image = fill_acceptance_image(image, center, offset, acceptance)

print ""

print "image.data"
print(image.data)
print "image.data[45:55, 45:55]"
print(image.data[45:55, 45:55])

print ""

#plot (filled) image
print "here comes the plot"
fig = plt.figure()
ax = fig.add_subplot(111)

ax.imshow(image.data[:,:], origin='lower')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('filled image')

plt.draw() #draw plot

print ""

#test: sum of the image (weigthed by the pixel areas) should be equal to the integral of the radial acceptance

#initialize WCS to the header of the image
w = WCS(image.header)

#define grids of pixel coorinates
nx, ny = image.data.shape
xpix_coord_grid = np.zeros(image.shape)
for y in range(0, ny):
    xpix_coord_grid[0:nx, y] = np.arange(0, nx)
print "xpix_coord_grid"
print(xpix_coord_grid) #debug
ypix_coord_grid = np.zeros(image.shape)
for x in range(0, nx):
    ypix_coord_grid[x, 0:nx] = np.arange(0, ny)
print "ypix_coord_grid"
print(ypix_coord_grid) #debug

print ""

#calculate pixel coordinates (in world coordinates)
coord = pixel_to_skycoord(xpix_coord_grid, ypix_coord_grid, w, 1)
print "coord"
print(coord) #debug

print ""

#pixel area = delta x * delta y * cos(zenith)
# zenith is either declination or latitude
pix_area_grid = np.cos(coord.l.to(u.radian))*x_pix_size.to(u.radian)*y_pix_size.to(u.radian)
print "pix_area_grid"
print(pix_area_grid) #debug

#sum image, weighted by pixel sizes (i.e. calculate integral of the image)
image_int_grid = image.data*pix_area_grid
image_int = image_int_grid.sum()
print "image_int:", image_int
print "image_int:", image_int.to(u.degree**2)

#integrate acceptance (i.e. gaussian function)
#1st integrate in r
acc_int = Quantity(quad(radial_gaussian_2D, Angle(0.*u.degree).to(u.radian).value,  Angle(10.*u.degree).to(u.radian).value, args=(sigma.to(u.radian).value))*u.radian)
print "acc_int:", acc_int
acc_int_value = acc_int[0]
print "acc_int_value:", acc_int_value
#2nd integrate in phi: phi range [0, 2pi)
acc_int_value *= Angle(2.*np.pi*u.radian)
print "acc_int_value:", acc_int_value
print "acc_int_value:", acc_int_value.to(u.degree**2)

#check sum of the image:
# int ~= sum (bin content * bin size)
# this is true if the pixelation is not too coarse (i.e. bin size small enough w.r.t. dimensions of the structure in the image (in this case the gaussian sigma))
epsilon = 1.e-4
assert abs(image_int.to(u.rad**2).value - acc_int_value.to(u.rad**2).value) < epsilon, "image integral not compatible with radial acceptance integral"

#TODO: save fits and image (eps, pdf, png)! (and check fits in ds9/fv)!!!!

plt.show() #don't quit at the end

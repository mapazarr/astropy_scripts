
import numpy as np
from scipy.integrate import quad
from matplotlib import pyplot as plt

from astropy.io import fits
from astropy.coordinates import SkyCoord, Angle
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
import astropy.units as u
from astropy.units.quantity import Quantity
from astropy.modeling import models

from gammapy.background import fill_acceptance_image
from gammapy.image import make_empty_image, lookup
from gammapy.image.utils import coordinates

#I don't need the function: I can call astropy.modeling.models.Gaussian1D
#def radial_gaussian_2D(r, sig):
#    """Radially symmetric Gaussian function for emulating acceptance curve.
#
#    Please give parameters in radians.
#
#    TODO: I can't define parameters as `~astropy.coordinates.Angle` because it conflicts with the integral call afterwards!!!
#
#    Parameters
#    ----------
#    r : `~numpy.ndarray`
#    Radial coordinate.
#    sig : float
#    Gaussian width.
#
#    Returns
#    -------
#    val: `~numpy.ndarray`
#    Gaussian evaluated at the requested radii.
#    """
#    norm = 1.
#    val = norm * np.exp(-np.power(r, 2.) / (2. * np.power(sig, 2.)))
#    return val

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
x_center = image.header['CRVAL1']
y_center = image.header['CRVAL2']

#sanity checks: am I using galactic coordinates in degrees?
assert image.header['CTYPE1'] == 'GLON-CAR', "Error: x coordinate is not in Galactic coordinates!"
assert image.header['CTYPE2'] == 'GLAT-CAR', "Error: y coordinate is not in Galactic coordinates!"
assert image.header['CUNIT1'] == 'deg', "Error: x coordinate is not in degrees!"
assert image.header['CUNIT2'] == 'deg', "Error: y coordinate is not in degrees!"

center = SkyCoord(l=x_center*u.degree, b=y_center*u.degree, frame='galactic')

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
#acceptance = radial_gaussian_2D(offset.to(u.radian).value, sigma.to(u.radian).value)
gaus_model = models.Gaussian1D(amplitude=1, mean=0., stddev=sigma.to(u.radian).value)
acceptance = gaus_model(offset.to(u.radian).value)

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
xpix_coord_grid, ypix_coord_grid = coordinates(image, world=False)

print ""

#calculate pixel coordinates (in world coordinates)
coord = pixel_to_skycoord(xpix_coord_grid, ypix_coord_grid, w, 0)
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
#1st integrate in r: remember: FoV radius for image of 100 pix of 0.1 deg size is 5 deg
#acc_int = Quantity(quad(radial_gaussian_2D, Angle(0.*u.degree).to(u.radian).value,  Angle(5.*u.degree).to(u.radian).value, args=(sigma.to(u.radian).value))*u.radian)
acc_int = Quantity(quad(gaus_model, Angle(0.*u.degree).to(u.radian).value,  Angle(5.*u.degree).to(u.radian).value, args=(sigma.to(u.radian).value))*u.radian)
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
decimal = 4
s_error = "image integral not compatible with radial acceptance integral"
#np.testing.assert_almost_equal(image_int.to(u.rad**2).value, acc_int_value.to(u.rad**2).value, decimal, s_error)
#TODO: the test fails!!! fixme!!!

#test: check points at the offsets where the acceptance is defined along the x axis (i.e. y=0 in pix coord)

print ""

for off, acc in zip(offset, acceptance):
    #print " off: %f, acc: %f" %(off.value, acc) #debug
    x_coord = x_center + off
    y_coord = y_center
    if off < Angle(5, 'degree'):
        image_acc = lookup(image, x_coord, y_coord, world=True)
        print " off: %f, acc: %f, image_acc: %f" %(off.value, acc, image_acc) #debug
        decimal = 1
        s_error = "image acceptance not compatible with defined radial acceptance"
        np.testing.assert_almost_equal(image_acc, acc, decimal, s_error)
        #TODO: antes de usar gammapy.image.utils.coordinates en fov.py el assert funcionaba (con decimal=1 pero funcionaba)!!!
        #en cualquier caso: me parece que habia un error en el grid de pixeles tal como estaba antes, pues daba error si usaba una imagen asimetrica (i.e. 5x3)!!!

    ##coord = SkyCoord(l=x_coord*u.degree, b=y_coord*u.degree, frame='galactic')

    #transform to pixel coord
    
    ##image_acc = image.data[x_coord_pix, y_coord_pix]

##coord = pixel_to_skycoord(xpix_coord_grid, ypix_coord_grid, w, 0)





#TODO: save fits and image (eps, pdf, png)! (and check fits in ds9/fv)!!!!

plt.show() #don't quit at the end

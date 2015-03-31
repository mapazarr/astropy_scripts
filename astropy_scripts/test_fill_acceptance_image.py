
from gammapy.background import fill_acceptance_image
from gammapy.image import make_empty_image

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
import astropy.units as u

import numpy as np
from scipy.integrate import quad
from matplotlib import pyplot as plt

#gaussian function (for emulating radial acceptane)
def radial_gaussian_2D(r, sig):
    norm = 1./(sig * np.sqrt(2. * np.pi))
    norm *= 2. #due to definition only for r>=0
    norm /= (2.*np.pi) #due to raial symmetry
    return  norm * np.exp(-np.power(r, 2.) / (2. * np.power(sig, 2.)))

#create empty image
image = make_empty_image()
#image = make_empty_image(1000,1000)
#image = make_empty_image(50, 50)
#image = make_empty_image(10, 10)
#image = make_empty_image(5, 10)

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
x_pix_size = abs(image.header['CDELT1'])
y_pix_size = abs(image.header['CDELT2'])

#sanity checks: am I using galactic coordinates in degrees?
assert image.header['CTYPE1'] == 'GLON-CAR', "Error: x coordinate is not in Galactic coordinates!"
assert image.header['CTYPE2'] == 'GLAT-CAR', "Error: y coordinate is not in Galactic coordinates!"
assert image.header['CUNIT1'] == 'deg', "Error: x coordinate is not in degrees!"
assert image.header['CUNIT2'] == 'deg', "Error: y coordinate is not in degrees!"

center = SkyCoord(l=x*u.degree, b=y*u.degree, frame='galactic')

print "center"
print(center)

print ""

#define radial acceptance and offset angles
##offset = np.arange(0., 3., 0.1)
offset = np.arange(0., 30., 0.1)
##offset = np.arange(0., 30., 0.01)
##offset = np.arange(0., 300., 0.1)

#acceptance = np.empty([0])
#sigma = 1.0 #gaussian width
#for off in offset :
#    acceptance = np.append(acceptance, [gaussian(off, sigma)])
acceptance = np.zeros_like(offset)
sigma = 1.0 #gaussian width
acceptance = radial_gaussian_2D(offset, sigma)

print "offset"
print(offset)
print "acceptance"
print(acceptance)

#fill acceptance in the image
#alternative: create a new object, instead of overwriting the empty image
image = fill_acceptance_image(image, center, offset, acceptance)

#TODO: check that the acceptance is defined over the whole image!!!
#      offset should go beyond map dimensions!!!

#check also: plot_background_model.py from Christoph!!!

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
print(xpix_coord_grid) #debug
ypix_coord_grid = np.zeros(image.shape)
for x in range(0, nx):
    ypix_coord_grid[x, 0:nx] = np.arange(0, ny)
print(ypix_coord_grid) #debug

#calculate pixel coordinates (in world coordinates)
coord = pixel_to_skycoord(xpix_coord_grid, ypix_coord_grid, w, 1)
print(coord) #debug

#pixel area = delta x * delta y * cos(zenith)
# zenith is either declination or latitude
pix_area_grid = np.cos(coord.l)*x_pix_size*y_pix_size
print(pix_area_grid) #debug

#sum image, weighted by pixel sizes (i.e. calculate integral of the image)
image_int_grid = image.data*pix_area_grid
image_int = image_int_grid.sum()
print "image_int:", image_int

#integrate acceptance (i.e. gaussian function)
#acc_int = quad(radial_gaussian_2D, 0., float("inf"), args=(sigma))
acc_int = quad(radial_gaussian_2D, 0., 10., args=(sigma))
print "acc_int:", acc_int
acc_int_value = acc_int[0]
print "acc_int_value:", acc_int_value
#make it radial integral: phi range [0, 2pi)
acc_int_value *= 2.*np.pi
print "acc_int_value:", acc_int_value

#check sum of the image:
# int ~= sum (bin content * bin size)
# this is true if the pixelation is not too coarse (i.e. bin size small enough w.r.t. dimensions of the structure in the image (in this case the gaussian sigma))
#assert image_int == acc_int_value
epsilon = 1.e-4
assert image_int > acc_int_value - epsilon and image_int < acc_int_value + epsilon, "image integral not compatible with radial acceptance integral"

#TODO: save fits and image (eps, pdf, png)! (and check fits in ds9/fv)!!!!

plt.show() #don't quit at the end

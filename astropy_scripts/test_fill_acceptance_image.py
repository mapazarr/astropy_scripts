
from gammapy.background import fill_acceptance_image
from gammapy.image import make_empty_image

from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

import numpy as np
from matplotlib import pyplot as plt

import sys

#gaussian function (for emulating radial acceptane)
def gaussian(x, mu, sig):
    norm = 1./(sig * np.sqrt(2. * np.pi))
    return  norm * np.exp(-np.power(x - mu, 2.) / (2. * np.power(sig, 2.)))

#create empty image
#image = make_empty_image()
image = make_empty_image(10, 10)

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
print "image.data.shape: ", image.data.shape

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
x = image.header['CRVAL1']
y = image.header['CRVAL2']

#sanity checks: am I using galactic coordinates in degrees?
if image.header['CTYPE1'] != 'GLON-CAR' :
    print "error in x coordinate system!!!"
    sys.exit()
elif image.header['CTYPE2'] != 'GLAT-CAR' :
    print "error in y coordinate system!!!"
    sys.exit()
elif image.header['CUNIT1'] != 'deg' : 
   print "error in x coordinate unit!!!"
   sys.exit()
elif image.header['CUNIT2'] != 'deg' :
    print "error in y coordinate unit!!!"
    sys.exit()

center = SkyCoord(l=x*u.degree, b=y*u.degree, frame='galactic')

print "center"
print(center)

print ""

#define radial acceptance and offset angles
##offset = np.arange(0., 3., 0.1)
offset = np.arange(0., 30., 0.1)

#acceptance = np.empty([0])
#sigma = 1.0 #gaussian width
#for off in offset :
#    acceptance = np.append(acceptance, [gaussian(off, 0, sigma)])
acceptance = np.zeros_like(offset)
sigma = 1.0 #gaussian width
acceptance = gaussian(offset, 0., sigma)

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

#plot (filled) image
print "here comes the plot"
fig = plt.figure()
ax = fig.add_subplot(111)

ax.imshow(image.data[:,:], origin='lower')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('filled image')

plt.draw() #draw plot

#TODO: save fits and image! (and check fits in ds9/fv)!!!!
#      also for test pyfits!!!!!

#TODO: repo with my test scripts?!!!
#      or save them in ANALISIS scripts???!!!
#      better: do repo and create a dir in ANALISIS scripts with a ref to repo!!!

plt.show() # don't quit at the end

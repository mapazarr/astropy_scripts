
from gammapy.image import make_empty_image

from astropy.io import fits

import numpy as np
from matplotlib import pyplot as plt

#create empty image
image = make_empty_image()
#image = make_empty_image(10, 10)

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

#add random data to the image
#image.data = np.random.random((128,128))
image.data = np.random.random(image.data.shape)
print "image.data.shape: ", image.data.shape

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

#NOTE: los arrays estan definidos al reves: array[y:x], primero la coord "y" y luego la "x" en vez de al reves!!!
#TODO: a test to reflect this!!! (esto me afecta el test de test_fill_acceptance_image.py!!!

#TODO: save fits and image (eps, pdf, png)! (and check fits in ds9/fv)!!!!

plt.show()) #don't quit at the end

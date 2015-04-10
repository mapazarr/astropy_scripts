
import numpy as np
from scipy.integrate import quad
from matplotlib import pyplot as plt

from astropy.io import fits
from astropy.coordinates import SkyCoord, Angle
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
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

def test_fill_acceptance_image():

    #create empty image
    #I need odd number of pixels if I want to have the center in its own pixel
    n_pix_x = 101
    n_pix_y = 101
    #n_pix_x = 11
    #n_pix_y = 11
    #n_pix_x = 11
    #n_pix_y = 7
    #n_pix_x = 3
    #n_pix_y = 3
    bin_size = Angle(0.1, 'degree')
    image = make_empty_image(n_pix_x, n_pix_y, bin_size.to(u.degree).value,
                             xref=0, yref=0, fill=0,
                             proj='CAR', coordsys='GAL',
                             xrefpix=None, yrefpix=None, dtype='float32')

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

    #estimate FoV radius (maximum offset) of the image:

    fov_x = (n_pix_x*bin_size)/2.
    fov_y = (n_pix_y*bin_size)/2.

    print "FoV radius (%f, %f) deg" %(fov_x.to(u.degree).value, fov_y.to(u.degree).value)

    #initialize WCS to the header of the image
    w = WCS(image.header)

    x_center_pix, y_center_pix = skycoord_to_pixel(center, w, 0)

    print "center pixel: (%f, %f)" %(x_center_pix, y_center_pix) #debug: is the center pixel on its own pixel?
    print "center pixel: (%i, %i)" %(x_center_pix, y_center_pix)

    #define pixel sizes
    x_pix_size = Angle(abs(image.header['CDELT1']), 'degree')
    y_pix_size = Angle(abs(image.header['CDELT2']), 'degree')

    print ""

    #define radial acceptance and offset angles
    #using bin_size for the offset step makes the test comparison easier
    #offset = Angle(np.arange(0., 30., 0.1), 'degree')
    offset = Angle(np.arange(0., 30., bin_size.to(u.degree).value), 'degree')
    acceptance = np.zeros_like(offset)
    sigma = Angle(1.0, 'degree') #gaussian width
    #acceptance = radial_gaussian_2D(offset.to(u.radian).value, sigma.to(u.radian).value)
    gaus_model = models.Gaussian1D(amplitude=1, mean=0., stddev=sigma.to(u.radian).value)
    acceptance = gaus_model(offset.to(u.radian).value)

    print "offset"
    print(offset)
    print "offset in deg"
    print offset.to(u.degree).value
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
    #1st integrate in r: remember: integrate only to up to image FoV
    #acc_int = Quantity(quad(radial_gaussian_2D, Angle(0., 'degree').to(u.radian).value,  Angle(5., 'degree').to(u.radian).value, args=(sigma.to(u.radian).value)), 'radian')
    acc_int = Quantity(quad(gaus_model, Angle(0., 'degree').to(u.radian).value,  Angle(5., 'degree').to(u.radian).value, args=(sigma.to(u.radian).value)), 'radian')
    print "acc_int:", acc_int
    acc_int_value = acc_int[0]
    print "acc_int_value:", acc_int_value
    #2nd integrate in phi: phi range [0, 2pi)
    acc_int_value *= Angle(2.*np.pi, 'radian')
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
        if off < fov_x:
            image_acc = lookup(image, x_coord, y_coord, world=True)
            print " off: %f, acc: %f, image_acc: %f" %(off.value, acc, image_acc) #debug
            decimal = 1
            #s_error = "image acceptance not compatible with defined radial acceptance"
            #np.testing.assert_almost_equal(image_acc, acc, decimal, s_error)
            #TODO: antes de usar gammapy.image.utils.coordinates en fov.py el assert funcionaba (con decimal=1 pero funcionaba)!!!
            #en cualquier caso: me parece que habia un error en el grid de pixeles tal como estaba antes, pues daba error si usaba una imagen asimetrica (i.e. 5x3)!!!
            #el problema era que los arrays estan definidos al reves: array[y:x], primero la coord "y" y luego la "x" en vez de al reves!!!
            #pero el test sigue sin funcionar!!! fixme!!!

    #test: check points at the offsets where the acceptance is defined along the x axis (i.e. y=0 in pix coord)
    #      another approach to get the data: follow te same functions to get pix coord as in the function to test

    # define grids of pixel coorinates
    xpix_coord_grid, ypix_coord_grid = coordinates(image, world=False)

    print ""

    # calculate pixel offset from center (in world coordinates)
    coord = pixel_to_skycoord(xpix_coord_grid, ypix_coord_grid, w, 0)
    pix_off = coord.separation(center)

    print "pix_off in deg"
    print (pix_off.to(u.degree).value)

    print "image.data"
    print (image.data)

    print ""

    #plot pixel offsets
    print "here comes the plot"
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.imshow(pix_off.to(u.degree).value[:,:], origin='lower')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('pixel offset in deg')

    plt.draw() #draw plot

    print ""

    # x axis (i.e. y=0 in pix coord) defined in the array positions [y_center_pix,:]
    # only interested in semi axis, so [y_center_pix, x_center_pix:]
    #WARNING: los arrays estan definidos al reves: array[y:x], primero la coord "y" y luego la "x" en vez de al reves!!!
    print "x_center_pix = ", x_center_pix #debug
    print "int(round(x_center_pix)) = ", int(round(x_center_pix)) #debug
    print "y_center_pix = ", y_center_pix #debug
    print "int(round(y_center_pix)) = ", int(round(y_center_pix)) #debug
    pix_off_x_axis = pix_off[int(round(y_center_pix)), int(round(x_center_pix)):]
    image.data_x_axis = image.data[int(round(y_center_pix)), int(round(x_center_pix)):]

    print "pix_off_x_axis in deg"
    print (pix_off_x_axis.to(u.degree).value)

    print "image.data_x_axis"
    print(image.data_x_axis)

    print ""

    # cut offset and acceptance arrays to match image size
    # this is only valid if the offset step matches the pixel size!!!
    n = pix_off_x_axis.size
    offset_cut = offset[0:n]
    acceptance_cut = acceptance[0:n]

    print "offset_cut in deg"
    print (offset_cut.to(u.degree).value)

    print "acceptance_cut"
    print(acceptance_cut)

    print ""

    decimal = 4
    s_error = "image acceptance not compatible with defined radial acceptance"
    np.testing.assert_almost_equal(image.data_x_axis, acceptance_cut, decimal, s_error)

    print ""

    #TODO: try not to transform offsets to deg and get value until the end!!!! (in case I need it... maybe i don't!!!)

    #TODO: offset angle limitation is hard coded!!!!

    #TODO: save fits and image (eps, pdf, png)! (and check fits in ds9/fv)!!!!

    print "end"

    plt.show() #don't quit at the end

if __name__=="__main__":
    test_fill_acceptance_image()

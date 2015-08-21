from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
from tempfile import NamedTemporaryFile
import numpy as np
from numpy.testing import assert_allclose
import matplotlib.pyplot as plt
from astropy.units import Quantity
from astropy.coordinates import Angle
from astropy.tests.helper import assert_quantity_allclose
from gammapy.background import Cube
from gammapy import datasets
from gammapy_bg_cube_models_utilities import CubeUtils
from gammapy.datasets.make import make_test_bg_cube_model

GRAPH_DEBUG = 0
USE_TEMP_FILES = 0
CACHE = 1 # set to 0 to ignore downloaded files in the cache
          # (or use rm ~/.astropy/cache)

def bg_cube_model_plots(filename):
    """Plot cube background model and store it in fits.

    The 'image' format file can be viewed with ds9.
    """
    bg_cube_model = Cube.read(filename, format='table', scheme='bg_cube')

    ##bg_cube_model.plot_images() # old code
    CubeUtils.plot_images(bg_cube_model)
    CubeUtils.plot_images(bg_cube_model, energy=Quantity(2., 'TeV'))
    #bg_cube_model.plot_image(energy=Quantity(1., 'TeV'))
    bg_cube_model.plot_image(energy=Quantity(2., 'TeV'))
    #bg_cube_model.plot_image(energy=Quantity(3., 'TeV'))
    #bg_cube_model.plot_image(energy=Quantity(1.53216636, 'TeV'))
    #bg_cube_model.plot_image(energy=Quantity(2.27553916, 'TeV'))
    #bg_cube_model.plot_image(energy=Quantity(80., 'TeV'))
    #bg_cube_model.plot_image(energy=Quantity([1.], 'TeV'))
    #bg_cube_model.plot_image(energy=Quantity([2.], 'TeV'))
    #bg_cube_model.plot_image(energy=Quantity([1., 2.], 'TeV'))

    ##bg_cube_model.plot_spectra() # old code
    CubeUtils.plot_spectra(bg_cube_model, format='mosaic')
    CubeUtils.plot_spectra(bg_cube_model, format='stack')
    CubeUtils.plot_spectra(bg_cube_model, coord=Angle([0., 0.], 'degree'))
    bg_cube_model.plot_spectrum(coord=Angle([0., 0.], 'degree'))
    #bg_cube_model.plot_spectrum(coord=Angle([0., 2.], 'degree'))
    #bg_cube_model.plot_spectrum(coord=Angle([-5., 0.], 'degree'))
    #bg_cube_model.plot_spectrum(coord=Angle([0., -5.], 'degree'))
    #bg_cube_model.plot_spectrum(coord=Angle([-5., -5.], 'degree'))
    #bg_cube_model.plot_spectrum(coord=Angle([[0., 0.]], 'degree'))
    #bg_cube_model.plot_spectrum(coord=Angle([[0., 0.], [1., 1.]], 'degree'))

    outfile = 'cube_background_model.fits'
    bg_cube_model.write(outfile, format='image', clobber=True) # overwrite


def test_cube_class(filename):
    """Testing the tests for gammapy.
    """
    # test shape and scheme of bg cube when reading a file
    scheme='bg_cube'
    bg_cube_model = Cube.read(filename, format='table', scheme=scheme)
    print("bg_cube_model.data.shape", bg_cube_model.data.shape)
    print("len(bg_cube_model.data.shape)", len(bg_cube_model.data.shape))
    assert len(bg_cube_model.data.shape) == 3
    assert bg_cube_model.data.shape == (len(bg_cube_model.energy_edges) - 1,
                                              len(bg_cube_model.coordy_edges) - 1,
                                              len(bg_cube_model.coordx_edges) - 1)
    assert bg_cube_model.scheme == scheme

    # example how to access data in cube
    energy_bin = bg_cube_model.find_energy_bin(energy=Quantity(2., 'TeV'))
    coord_bin = bg_cube_model.find_coord_bin(coord=Angle([0., 0.], 'degree'))
    bg_cube_model.data[energy_bin, coord_bin[1], coord_bin[0]]


    # test image plot:
    # test bg rate values plotted for image plot of energy bin conaining E = 2 TeV
    energy = Quantity(2., 'TeV')
    ax_im = bg_cube_model.plot_image(energy)
    plt.draw()
    if GRAPH_DEBUG:
        plt.show() # wait until image is closed
    # get plot data (stored in the image)
    image_im = ax_im.get_images()[0]
    plot_data = image_im.get_array()
    # I need also the bin ids!!!
    plot_shape = image_im.get_array().shape
    plot_extent = image_im.get_extent()
    print("plot data")
    print(plot_data)
    print("plot shape", plot_shape)
    print("plot extent", plot_extent)

    # get data from bg model object to compare
    energy_bin = bg_cube_model.find_energy_bin(energy)
    model_data = bg_cube_model.data[energy_bin]
    # TODO: get also coord (x,y coord) of the bins!!!
    print("model data")
    print(model_data)

    # test if both arrays are equal
    assert_allclose(plot_data, model_data.value)

    # test spectrum plot:
    # test bg rate values plotted for spectrum plot of detector bin conaining coord (0, 0) deg (center)
    coord = Angle([0., 0.], 'degree')
    #coord = Angle([0., 2.], 'degree')
    ax_spec = bg_cube_model.plot_spectrum(coord)
    plt.draw()
    if GRAPH_DEBUG:
        plt.show() # wait until image is closed
    # get plot data (stored in the line)
    plot_data = ax_spec.get_lines()[0].get_xydata()
    print("plot data")
    print(plot_data)

    # get data from bg model object to compare
    coord_bin = bg_cube_model.find_coord_bin(coord)
    model_data = bg_cube_model.data[:, coord_bin[1], coord_bin[0]]
    # TODO: get also energies (x coord) of the points!!!
    print("model data")
    print(model_data)

    # test if both arrays are equal
    assert_allclose(plot_data[:,1], model_data.value)

    # test write (save)
    # test if values are correct in the saved file: compare both files
    bg_cube_model_1 = bg_cube_model
    if not USE_TEMP_FILES:
        outfile = 'bg_cube_model.fits'
    else:
        outfile = NamedTemporaryFile(suffix='.fits').name
    print("Writing file {}".format(outfile))
    bg_cube_model_1.write(outfile, format='table', clobber=True) # overwrite
    bg_cube_model_2 = Cube.read(outfile, format='table', scheme='bg_cube')
    assert_quantity_allclose(bg_cube_model_2.data, bg_cube_model_1.data)
    assert_quantity_allclose(bg_cube_model_2.coordx_edges, bg_cube_model_1.coordx_edges)
    assert_quantity_allclose(bg_cube_model_2.coordy_edges, bg_cube_model_1.coordy_edges)
    assert_quantity_allclose(bg_cube_model_2.energy_edges, bg_cube_model_1.energy_edges)

    # test read/write image file
    if not USE_TEMP_FILES:
        outfile = 'bg_cube_model_image.fits'
    else:
        outfile = NamedTemporaryFile(suffix='.fits').name
    print("Writing file {}".format(outfile))
    bg_cube_model.write(outfile, format='image', clobber=True) # overwrite
    # test: this function should produce:
    #  - a file exactly as the other method
    #  - a cube exactly as the other method (once the file is read with the read function)

    # test if values are correct in the saved file: compare both files
    bg_cube_model_1 = bg_cube_model
    bg_cube_model_2 = Cube.read(outfile, format='image', scheme='bg_cube')
    assert_quantity_allclose(bg_cube_model_2.data, bg_cube_model_1.data)
    assert_quantity_allclose(bg_cube_model_2.coordx_edges, bg_cube_model_1.coordx_edges)
    assert_quantity_allclose(bg_cube_model_2.coordy_edges, bg_cube_model_1.coordy_edges)
    assert_quantity_allclose(bg_cube_model_2.energy_edges, bg_cube_model_1.energy_edges)

    # TODO: test cubes/plots with asymmetric shape (x_bins != y_bins) !!!

    # TODO: doc for plots, example plots (revert read/write changes? -> no) !!!
    #       my scripts: plot stack + clean code
    #       review "use_plot_something.py"

    # TODO: clean up after tests (remove created files) !!!!!!!!!!!!!!!!!!!!!! (add option/flag to aventually keep them)

    plt.show() #don't quit at the end


def test_make_test_bg_cube_model(debug=False):
    """
    Testing make_test_bg_cube_model.
    """
    # call the function without arguments
    bg_cube_model = make_test_bg_cube_model()

    # plots (images and spectra)
    if debug:
        #make plots
        CubeUtils.plot_images(bg_cube_model)
        CubeUtils.plot_spectra(bg_cube_model, format='stack') # slow!

    outfile = 'test_bg_cube_model_table.fits'
    print("Writing file {}".format(outfile))
    bg_cube_model.write(outfile, format='table', clobber=True) # overwrite

    outfile = 'test_bg_cube_model_image.fits'
    print("Writing file {}".format(outfile))
    bg_cube_model.write(outfile, format='image', clobber=True) # overwrite

    # make a cube bg model with non-equal axes
    ndetx_bins = 1
    ndety_bins = 2
    nenergy_bins = 3
    bg_cube_model = make_test_bg_cube_model(ndetx_bins=ndetx_bins, ndety_bins=ndety_bins, nenergy_bins=nenergy_bins)

    # plots (images and spectra)
    if debug:
        #make plots
        CubeUtils.plot_images(bg_cube_model)
        CubeUtils.plot_spectra(bg_cube_model, format='stack') # slow!

    outfile = 'test_bg_cube_model_xydiff_table.fits'
    print("Writing file {}".format(outfile))
    bg_cube_model.write(outfile, format='table', clobber=True) # overwrite

    outfile = 'test_bg_cube_model_xydiff_image.fits'
    print("Writing file {}".format(outfile))
    bg_cube_model.write(outfile, format='image', clobber=True) # over

    # test shape of cube bg model
    assert len(bg_cube_model.data.shape) == 3
    assert bg_cube_model.data.shape == (nenergy_bins, ndety_bins, ndetx_bins)

    # make masked bg model
    bg_cube_model = make_test_bg_cube_model(apply_mask=True)

    # plots (images and spectra)
    if debug:
        #make plots
        CubeUtils.plot_images(bg_cube_model)
        #CubeUtils.plot_spectra(bg_cube_model, format='stack') # slow! #this should fail because of 0 plots in log axis!

    outfile = 'test_bg_cube_model_masked_table.fits'
    print("Writing file {}".format(outfile))
    bg_cube_model.write(outfile, format='table', clobber=True) # overwrite

    outfile = 'test_bg_cube_model_masked_image.fits'
    print("Writing file {}".format(outfile))
    bg_cube_model.write(outfile, format='image', clobber=True) # overwrite

    # test that values with (x, y) > (0, 0) are zero
    x_points = Angle(np.arange(5), 'degree') + Angle(0.01, 'degree')
    y_points = Angle(np.arange(5), 'degree') + Angle(0.01, 'degree')
    e_points = bg_cube_model.energy_bin_centers
    x_points, y_points, e_points = np.meshgrid(x_points, y_points, e_points,
                                               indexing='ij')
    det_bin_index = bg_cube_model.find_coord_bin(Angle([x_points, y_points]))
    e_bin_index = bg_cube_model.find_energy_bin(e_points)
    bg = bg_cube_model.data[e_bin_index, det_bin_index[1], det_bin_index[0]]

    #assert that values are 0
    assert_quantity_allclose(bg, Quantity(0., bg.unit))

def test_remote_data():
    """Testing the use of remote data in gammapy.
    """

    # test local path
    # checks (localy) in gammapy/gammapy/datasets/data
    print()
    local_path = datasets.get_path('fermi/fermi_counts.fits.gz')
    print("local_path", local_path)

    # test remote path
    # checks (remotely) in gammapy/gammapy-extra/datasets
    print()
    remote_path = datasets.get_path('vela_region/counts_vela.fits', location='remote', cache=CACHE)
    print("remote_path", remote_path)

    # test local path
    print()
    # checks (localy) in gammapy/gammapy/datasets/data
    #local_path = datasets.get_path('../../background/tests/data/bg_test.fits')
    local_path = datasets.get_path('hess/run_0023037_hard_eventlist.fits.gz')
    print("local_path", local_path)

    # test remote path
    print()
    # checks (remotely) in gammapy/gammapy-extra/datasets
    #remote_path = datasets.get_path('bg_cube_model_test.fits', location='remote', cache=CACHE)
    remote_path = datasets.get_path('../test_datasets/background/bg_cube_model_test1.fits', location='remote', cache=CACHE)
    print("remote_path", remote_path)


if __name__ == '__main__':

    # create a list of files to test
    filenames = []

    ##DIR = '/home/mapaz/HESS/fits_data/pa_fits_prod01/pa/Model_Deconvoluted_Prod26/Mpp_Std/background/
    ##DIR = '/home/mapaz/HESS/fits_data/pa_fits_prod02/pa/Model_Deconvoluted_Prod26/Mpp_Std/background/
    #DIR = '/home/mapaz/astropy/testing_cube_bg_michael_mayer/pa_fits_prod01/pa_fits_prod01_results/background_smooth_models/'
    DIR = '/home/mapaz/astropy/testing_cube_bg_michael_mayer/pa_fits_prod02/background/'
    filename = DIR + 'hist_alt3_az0.fits.gz'
    filenames.append(filename)

    #DIR = '/home/mapaz/astropy/development_code/gammapy/gammapy/background/tests/data/'
    #filename = DIR + 'bg_test.fits'
    #filename = DIR + 'bkgcube.fits' # not supported!
    filename = '../test_datasets/background/bg_cube_model_test1.fits'
    filename = datasets.get_path(filename, location='remote', cache=CACHE)
    #DIR = '/home/mapaz/astropy/development_code/gammapy-extra/test_datasets/background/'
    #filename = DIR + 'bg_cube_model_test.fits'
    filenames.append(filename)

    DIR = '/home/mapaz/astropy/development_code/astropy_scripts/astropy_scripts/'
    filename = DIR + 'test_bg_cube_model.fits'
    #filenames.append(filename)

    for filename in filenames:
        print()
        print("filename: {}".format(filename))

        # call tests
        #bg_cube_model_plots(filename) # takes long! (many plots/files created!)
        test_cube_class(filename)

    # call tests
    test_make_test_bg_cube_model(debug=False)
    #test_make_test_bg_cube_model(debug=True)
    #test_remote_data()

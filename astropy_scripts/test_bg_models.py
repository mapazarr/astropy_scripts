from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
from astropy.units import Quantity
from astropy.coordinates import Angle
from gammapy.background.models import CubeBackgroundModel

def plot_example():
    """Plot background model and store as cube so that it can viewed with ds9.
    """
    #DIR = '/Users/deil/work/_Data/hess/HESSFITS/pa/Model_Deconvoluted_Prod26/Mpp_Std/background/'
    DIR = '/home/mapaz/astropy/testing_cube_bg_michael_mayer/background/'
    filename = DIR + 'hist_alt3_az0.fits.gz'
    bg_model = CubeBackgroundModel.read(filename)

    bg_model.plot_images()
    #bg_model.plot_images(energy=Quantity(1., 'TeV'))
    #bg_model.plot_images(energy=Quantity(2., 'TeV'))
    #bg_model.plot_images(energy=Quantity(3., 'TeV'))
    #bg_model.plot_images(energy=Quantity(1.53216636, 'TeV'))
    #bg_model.plot_images(energy=Quantity(2.27553916, 'TeV'))
    #bg_model.plot_images(energy=Quantity(80., 'TeV'))
    #bg_model.plot_images(energy=Quantity([1.], 'TeV'))
    #bg_model.plot_images(energy=Quantity([2.], 'TeV'))
    #bg_model.plot_images(energy=Quantity([1., 2.], 'TeV'))

    bg_model.plot_spectra()
    #bg_model.plot_spectra(det=Angle([0., 0.], 'degree'))
    #bg_model.plot_spectra(det=Angle([0., 2.], 'degree'))
    #bg_model.plot_spectra(det=Angle([-5., 0.], 'degree'))
    #bg_model.plot_spectra(det=Angle([0., -5.], 'degree'))
    #bg_model.plot_spectra(det=Angle([-5., -5.], 'degree'))
    #bg_model.plot_spectra(det=Angle([[0., 0.]], 'degree'))
    #bg_model.plot_spectra(det=Angle([[0., 0.], [1., 1.]], 'degree'))

    bg_model.write_cube('cube_background_model.fits')

def gammapy_tests():
    """Testing the tests for gammapy.
    """
    # test shape of bg cube when reading a file
    DIR = '/home/mapaz/astropy/development_code/gammapy/gammapy/background/tests/data/'
    #filenames = ['bg_test.fits', 'bkgcube.fits']
    filenames = ['bg_test.fits']
    for filename in filenames:
        print()
        print("filename: {}".format(filename))
        filename = DIR + filename
        bg_cube_file = CubeBackgroundModel.read(filename)
        print("bg_cube_file.background.shape", bg_cube_file.background.shape)
        print("len(bg_cube_file.background.shape)", len(bg_cube_file.background.shape))
        assert len(bg_cube_file.background.shape) == 3
        #import IPython; IPython.embed()

if __name__ == '__main__':
    #plot_example()
    gammapy_tests()

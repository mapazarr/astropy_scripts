from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
import matplotlib.pyplot as plt
from astropy.units import Quantity
from astropy.coordinates import Angle
from gammapy.background.models import CubeBackgroundModel

def plot_example():
    """Plot background model and store as cube so that it can viewed with ds9.
    """
    #DIR = '/Users/deil/work/_Data/hess/HESSFITS/pa/Model_Deconvoluted_Prod26/Mpp_Std/background/'
    #DIR = '/home/mapaz/astropy/testing_cube_bg_michael_mayer/background/'
    #filename = DIR + 'hist_alt3_az0.fits.gz'
    DIR = '/home/mapaz/astropy/development_code/gammapy/gammapy/background/tests/data/'
    filename = DIR + 'bg_test.fits'
    bg_model = CubeBackgroundModel.read(filename)

    bg_model.plot_images()
    #bg_model.plot_images(energy=Quantity(1., 'TeV'))
    bg_model.plot_images(energy=Quantity(2., 'TeV'))
    #bg_model.plot_images(energy=Quantity(3., 'TeV'))
    #bg_model.plot_images(energy=Quantity(1.53216636, 'TeV'))
    #bg_model.plot_images(energy=Quantity(2.27553916, 'TeV'))
    #bg_model.plot_images(energy=Quantity(80., 'TeV'))
    #bg_model.plot_images(energy=Quantity([1.], 'TeV'))
    #bg_model.plot_images(energy=Quantity([2.], 'TeV'))
    #bg_model.plot_images(energy=Quantity([1., 2.], 'TeV'))

    bg_model.plot_spectra()
    bg_model.plot_spectra(det=Angle([0., 0.], 'degree'))
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
        bg_cube_model = CubeBackgroundModel.read(filename)
        print("bg_cube_model.background.shape", bg_cube_model.background.shape)
        print("len(bg_cube_model.background.shape)", len(bg_cube_model.background.shape))
        assert len(bg_cube_model.background.shape) == 3
        #import IPython; IPython.embed()

    # test image plot:
    # test bg rate values plotted for image plot of energy bin conaining E = 2 TeV
    fig_image, ax_im, image_im = bg_cube_model.plot_images(energy=Quantity(2., 'TeV'))
    plt.draw()#TODO: use graph debug!!!! (or not?)
    print("plot data")
    print(image_im.get_array())
    # I need also the bin ids!!!
    print("plot shape", image_im.get_array().shape)
    print("plot extent", image_im.get_extent())




    # test spectrum plot:
    # test bg rate values plotted for spectrum plot of detector bin conaining det (0, 0) deg (center)
    fig_spec, ax_spec, image_spec = bg_cube_model.plot_spectra(det=Angle([0., 0.], 'degree'))
    plt.draw()#TODO: use graph debug!!!! (or not?)
    print("plot data")
    #print(ax.get_lines()[0].get_xydata())
    print(ax_spec.get_lines()[0].get_xydata())



    plt.show()#TODO: use graph debug!!!! (or not?)

    #TODO: creo que el archivo de ejemplo tiene detx dety en radianes!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!(y yo con las unidades hardcoded :-D)!!!!
    # mirar si ese archivo tiene las unidades por algun sitio....!!!!!!!!!!!!!!!!!!!

if __name__ == '__main__':
    #plot_example()
    gammapy_tests()

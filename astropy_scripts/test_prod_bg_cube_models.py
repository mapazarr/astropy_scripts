from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
import os
import matplotlib.pyplot as plt
from astropy.units import Quantity
from astropy.coordinates import Angle
from gammapy.scripts import make_bg_cube_models
from gammapy.background import Cube
from gammapy.obs import ObservationGroups


DEBUG = 2 # 0: no output, 1: output, 2: NOTHING, 3: more verbose
GRAPH_DEBUG = 1 # 0: no plots, 1: make plots, 2: wait between steps (bins), 3: draw 3D scatter plots (not implemented)
CLEAN_WORKING_DIR = 1 # remove existing observation and bg cube model files

HESSFITS_MPP = '/home/mapaz/astropy/gammapy_tutorial/HESS_fits_data/pa/Model_Deconvoluted_Prod26/Mpp_Std'

def bg_models_debug_plots():

    # TODO: call plot_bg_cube_model_comparison !!!

    # read observation grouping
    indir = os.environ['PWD'] + '/'
    infile = indir + 'bg_observation_groups.ecsv'
    observation_groups = ObservationGroups.read('bg_observation_groups.ecsv')

    # loop over observation groups
    groups = observation_groups.list_of_groups
    if DEBUG:
        print()
        print("list of groups", groups)

    for group in groups:
        if DEBUG:
            print()
            print("group", group)

        # read bg cube model from file
        indir = os.environ['PWD'] + '/bg_cube_models/'
        infile = indir +\
                 'bg_cube_model_group{}_table.fits.gz'.format(group)
        # skip bins with no bg cube model file
        if not os.path.isfile(infile):
            print("WARNING, file not found: {}".format(infile))
            continue # skip the rest
        bg_cube_model = Cube.read(infile, format='table')

        fig, axes = plt.subplots(nrows=1, ncols=3)
        fig.set_size_inches(30., 8., forward=True)
        #plt.suptitle('altitude bin {0} azimuth bin {1}'.format(i_alt, i_az))
        plt.suptitle('group {}'.format(group))
        # TODO: it would be nice to get a nice string from an obs group!!!
        #       and in this case pack it in the figure title
        #       I think this also applies to the script for comparing plots!!!

        # plot images
        #bg_cube_model.plot_image(energy=Quantity(0.5, 'TeV'), ax=axes[0])
        bg_cube_model.plot_image(energy=Quantity(5., 'TeV'), ax=axes[0])
        bg_cube_model.plot_image(energy=Quantity(50., 'TeV'), ax=axes[1])

        # plot spectra
        bg_cube_model.plot_spectrum(det=Angle([0., 0.], 'degree'), ax=axes[2],
                                    style_kwargs=dict(color='blue',
                                                      label='(0, 0) deg'))
        bg_cube_model.plot_spectrum(det=Angle([2., 2.], 'degree'), ax=axes[2],
                                    style_kwargs=dict(color='red',
                                                      label='(2, 2) deg'))
        axes[2].set_title('')
        axes[2].legend()

        #plt.tight_layout()
        plt.draw()

        if GRAPH_DEBUG > 1:
            plt.show() # wait until image is closed

    if GRAPH_DEBUG:
        plt.show() # don't leave at the end


def test_make_bg_cube_models():
    """
    gammapy-background-cube -h
    gammapy-make_bg_cube_models -h
    gammapy-make_bg_cube_models -h
    gammapy-make_bg_cube_models /home/mapaz/astropy/gammapy_tutorial/HESS_fits_data/pa/Model_Deconvoluted_Prod26/Mpp_Std
    gammapy-make_bg_cube_models /home/mapaz/astropy/gammapy_tutorial/HESS_fits_data/pa/Model_Deconvoluted_Prod26/Mpp_Std --test True
    """
    # Need to make sure the working dir is clean, otherwise old
    # files could be mixed up in the new models!
    if CLEAN_WORKING_DIR:
        print("Cleaning working dir.")
        #command = "rm bg_observation_table* bg_cube_model_alt* -fr"
        #command = "rm bg_observation_table.fits.gz splitted_obs_list/ bg_cube_models/ -fr"
        command = "rm bg_observation_table.fits.gz bg_observation_groups.ecsv bg_observation_table_grouped.fits.gz bg_cube_models/ -fr"
        print(command)
        os.system(command)

    test = False
    if DEBUG > 1:
        # run fast (test mode)
        test = True

    make_bg_cube_models(fitspath=HESSFITS_MPP, test=test)

    if GRAPH_DEBUG:
        # check model: do some plots
        bg_models_debug_plots()


if __name__ == '__main__':
    test_make_bg_cube_models()

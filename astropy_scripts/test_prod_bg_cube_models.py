from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.units import Quantity
from astropy.coordinates import Angle
from gammapy.scripts import make_bg_cube_models
from gammapy.background import Cube
from gammapy.obs import ObservationGroups
from gammapy.datasets import make_test_dataset


DEBUG = 0 # 0: normal, 1: run fast (test mode)
GRAPH_DEBUG = 1 # 0: no plots, 1: make plots, 2: wait between steps (bins), 3: draw 3D scatter plots (not implemented)
CLEAN_WORKING_DIR = 1 # remove existing observation and bg cube model files
USE_DUMMY_DATA = 0 # to use dummy dataset

HESSFITS_MPP = '/home/mapaz/astropy/gammapy_tutorial/HESS_fits_data/pa/Model_Deconvoluted_Prod26/Mpp_Std'
DUMMYFITS = '/home/mapaz/astropy/development_code/astropy_scripts/astropy_scripts/' + 'test_dataset'
SCHEME = 'HESS'
OBSERVATORY_NAME = 'HESS' # in case USE_DUMMY_DATA is activated

def bg_cube_models_debug_plots(indir):
    """Make some debug plots of the generated background cube models.

    Parameters
    ----------
    indir : str
        Dir path where results are stored.
    """
    # TODO: call plot_bg_cube_model_comparison !!!

    print()
    print("#######################################")
    print("# Starting bg_cube_models_debug_plots #")
    print("#######################################")

    # read observation grouping
    infile = indir + '/bg_observation_groups.ecsv'
    observation_groups = ObservationGroups.read(infile)

    # loop over observation groups
    groups = observation_groups.list_of_groups
    print()
    print("list of groups", groups)

    for group in groups:
        print()
        print("group", group)

        # read bg cube model from file
        infile = indir +\
                 '/bg_cube_model_group{}_table.fits.gz'.format(group)
        # skip bins with no bg cube model file
        if not os.path.isfile(infile):
            print("WARNING, file not found: {}".format(infile))
            continue # skip the rest
        bg_cube_model = Cube.read(infile, format='table', scheme='bg_cube')

        fig, axes = plt.subplots(nrows=1, ncols=3)
        fig.set_size_inches(30., 8., forward=True)
        #plt.suptitle('altitude bin {0} azimuth bin {1}'.format(i_alt, i_az))
        #plt.suptitle('group {}'.format(group))
        group_info = observation_groups.info_group(group)
        plt.suptitle(group_info)
        # TODO: it would be nice to get a nice string from an obs group!!!
        #       and in this case pack it in the figure title
        #       I think this also applies to the script for comparing plots!!!

        # plot images
        #bg_cube_model.plot_image(energy=Quantity(0.5, 'TeV'), ax=axes[0])
        bg_cube_model.plot_image(energy=Quantity(5., 'TeV'), ax=axes[0])
        bg_cube_model.plot_image(energy=Quantity(50., 'TeV'), ax=axes[1])

        # plot spectra
        bg_cube_model.plot_spectrum(coord=Angle([0., 0.], 'degree'), ax=axes[2],
                                    style_kwargs=dict(color='blue',
                                                      label='(0, 0) deg'))
        bg_cube_model.plot_spectrum(coord=Angle([2., 2.], 'degree'), ax=axes[2],
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
    gammapy-make-bg-cube-models -h
    gammapy-make-bg-cube-models /home/mapaz/astropy/gammapy_tutorial/HESS_fits_data/pa/Model_Deconvoluted_Prod26/Mpp_Std HESS bg_cube_models
    gammapy-make-bg-cube-models /home/mapaz/astropy/gammapy_tutorial/HESS_fits_data/pa/Model_Deconvoluted_Prod26/Mpp_Std HESS bg_cube_models --test
    gammapy-make-bg-cube-models /home/mapaz/astropy/gammapy_tutorial/HESS_fits_data/pa/Model_Deconvoluted_Prod26/Mpp_Std HESS bg_cube_models --test --overwrite
    gammapy-make-bg-cube-models /home/mapaz/astropy/gammapy_tutorial/HESS_fits_data/pa/Model_Deconvoluted_Prod26/Mpp_Std HESS bg_cube_models --a-la-michi
    """
    # Need to make sure the working dir is clean, otherwise old
    # files could be mixed up in the new models!
    if CLEAN_WORKING_DIR:
        print("Cleaning working dir.")
        #command = "rm bg_observation_table* bg_cube_model_alt* -fr"
        #command = "rm bg_observation_table.fits.gz splitted_obs_list/ bg_cube_models/ -fr"
        #command = "rm bg_observation_table.fits.gz bg_observation_groups.ecsv bg_observation_table_grouped.fits.gz bg_cube_models/ -fr"
        command = "rm bg_cube_models/ -fr"
        print(command)
        os.system(command)
        if USE_DUMMY_DATA:
          print("Deleting dummy data dir.")
          command = "rm test_dataset -fr"
          print(command)
          os.system(command)

    test = False
    if DEBUG:
        # run fast (test mode)
        test = True

    fits_path = HESSFITS_MPP
    #outdir = os.environ['PWD'] + '/bg_cube_models/'
    outdir = 'bg_cube_models'
    overwrite = False
    #a_la_michi = False
    a_la_michi = True

    if USE_DUMMY_DATA:
        # update fits path and generate dataset
        fits_path = DUMMYFITS
        n_obs = 10
        if DEBUG:
            # run fast (test mode)
            n_obs = 2
        datestart = None
        dateend = None
        #random_state = 'random-seed'
        random_state = np.random.RandomState(seed=0)
        #random_state = 0 # this is equivalent
        make_test_dataset(fits_path=fits_path,
                          overwrite=overwrite,
                          observatory_name=OBSERVATORY_NAME,
                          n_obs=n_obs,
                          datestart=datestart,
                          dateend=dateend,
                          random_state=random_state)

    make_bg_cube_models(fitspath=fits_path, scheme=SCHEME, outdir=outdir, overwrite=overwrite, test=test, a_la_michi=a_la_michi)

    if GRAPH_DEBUG:
        # check model: do some plots
        bg_cube_models_debug_plots(indir=outdir)


if __name__ == '__main__':
    test_make_bg_cube_models()

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.units import Quantity
from astropy.coordinates import Angle
from gammapy.background import Cube

GRAPH_DEBUG = 0
SAVE = 0

#input_dir1 = '/home/mapaz/astropy/working_dir/hess-host_scripts/pa_fits_prod01_results/bg_cube_models'
input_dir1 = '/home/mapaz/astropy/working_dir/hess-host_scripts/pa_fits_prod02_results/bg_cube_models'
#input_dir2 = '/home/mapaz/HESS/fits_data/pa_fits_prod01/pa/Model_Deconvoluted_Prod26/Mpp_Std/background'
input_dir2 = '/home/mapaz/HESS/fits_data/pa_fits_prod02/pa/Model_Deconvoluted_Prod26/Mpp_Std/background'

# TODO: addapt to new mods of background cube classes!!!
#       this tool could benefit from an interpolator method to find the correct obs group to compare!!! (TODO: elaborate more this TODO idea) !!!

# alt az bin IDs for comparison
alt_bin_ids = [7, 10, 13]
az_bin_ids = [0, 1]

# alt az bin edges definitions
altitude_edges = Angle([0, 20, 23, 27, 30, 33, 37, 40, 44, 49, 53, 58, 64, 72, 90], 'degree')
azimuth_edges = Angle([-90, 90, 270], 'degree')

# loop over alt-az bins
for i_alt in alt_bin_ids:
    print()
    alt_bin_edges = Angle([altitude_edges[i_alt], altitude_edges[i_alt + 1]])
    print("bin alt {}, edges = {}".format(i_alt, alt_bin_edges))

    for i_az in az_bin_ids:
        print()
        az_bin_edges = Angle([azimuth_edges[i_az], azimuth_edges[i_az + 1]])
        print("bin az {}, edges = {}".format(i_az, az_bin_edges))

        # get cubes
        filename1 = input_dir1 + '/bg_cube_model_alt' + str(i_alt) +\
                    '_az' + str(i_az) + '_table.fits'
        filename2 = input_dir2 + '/hist_alt' + str(i_alt) +\
                    '_az' + str(i_az) + '.fits.gz'
        print('filename1', filename1)
        print('filename2', filename2)
        bg_cube_model1 = Cube.read(filename1, format='table', scheme='bg_cube')
        bg_cube_model2 = Cube.read(filename2, format='table', scheme='bg_cube')

        # compare binning
        print("energy edges 1", bg_cube_model1.energy_edges)
        print("energy edges 2", bg_cube_model2.energy_edges)
        print("detector edges 1 Y", bg_cube_model1.coordy_edges)
        print("detector edges 2 Y", bg_cube_model2.coordy_edges)
        print("detector edges 1 X", bg_cube_model1.coordx_edges)
        print("detector edges 2 X", bg_cube_model2.coordx_edges)

        # plot
        fig, axes = plt.subplots(nrows=2, ncols=3)
        fig.set_size_inches(30., 15., forward=True)
        plt.suptitle('altitude bin {0}: {2}    azimuth bin {1}: {3}'.format(i_alt, i_az, alt_bin_edges, az_bin_edges))

        # plot images
        #  rows: similar energy bin
        #  cols: same file
        #bg_cube_model1.plot_image(energy=Quantity(0.5, 'TeV'), ax=axes[0, 0])
        bg_cube_model1.plot_image(energy=Quantity(5., 'TeV'), ax=axes[0, 0])
        axes[0, 0].set_title("model 1: {}".format(axes[0, 0].get_title()))
        bg_cube_model1.plot_image(energy=Quantity(50., 'TeV'), ax=axes[1, 0])
        axes[1, 0].set_title("model 1: {}".format(axes[1, 0].get_title()))
        #bg_cube_model2.plot_image(energy=Quantity(0.5, 'TeV'), ax=axes[0, 1])
        bg_cube_model2.plot_image(energy=Quantity(5., 'TeV'), ax=axes[0, 1])
        axes[0, 1].set_title("model 2: {}".format(axes[0, 1].get_title()))
        bg_cube_model2.plot_image(energy=Quantity(50., 'TeV'), ax=axes[1, 1])
        axes[1, 1].set_title("model 2: {}".format(axes[1, 1].get_title()))

        # plot spectra
        #  rows: similar det bin
        #  cols: compare both files
        bg_cube_model1.plot_spectrum(coord=Angle([0., 0.], 'degree'), ax=axes[0, 2],
                                     style_kwargs=dict(color='blue',
                                                       label='model 1'))
        spec_title1 = axes[0, 2].get_title()
        bg_cube_model2.plot_spectrum(coord=Angle([0., 0.], 'degree'), ax=axes[0, 2],
                                    style_kwargs=dict(color='red',
                                                      label='model 2'))
        spec_title2 = axes[0, 2].get_title()
        if spec_title1 != spec_title2:
            raise ValueError("Expected same det binning, but got \"{0}\" and \"{1}\"".format(spec_title1, spec_title2))
        else:
            axes[0, 2].set_title(spec_title1)
        axes[0, 2].legend()

        bg_cube_model1.plot_spectrum(coord=Angle([2., 2.], 'degree'), ax=axes[1, 2],
                                    style_kwargs=dict(color='blue',
                                                      label='model 1'))
        spec_title1 = axes[1, 2].get_title()
        bg_cube_model2.plot_spectrum(coord=Angle([2., 2.], 'degree'), ax=axes[1, 2],
                                    style_kwargs=dict(color='red',
                                                      label='model 2'))
        spec_title2 = axes[1, 2].get_title()
        if spec_title1 != spec_title2:
            raise ValueError("Expected same det binning, but got \"{0}\" and \"{1}\"".format(spec_title1, spec_title2))
        else:
            axes[1, 2].set_title(spec_title1)
        axes[1, 2].legend()

        if GRAPH_DEBUG:
            plt.show() # wait until image is closed

        if SAVE:
            outfile = "bg_cube_model_comparison_alt{0}_az{1}.png".format(i_alt, i_az)
            print('Writing {}'.format(outfile))
            fig.savefig(outfile)

plt.show() # don't leave at the end

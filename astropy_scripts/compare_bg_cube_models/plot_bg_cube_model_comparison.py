"""
Script to produce plots comparing 2 sets of background cube models.

Details in stringdoc of the plot_bg_cube_model_comparison function.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm
from astropy.units import Quantity
from astropy.coordinates import Angle
from astropy.table import Table
from astropy.io import ascii
from gammapy.background import CubeBackgroundModel
from gammapy.obs import ObservationGroups, ObservationGroupAxis

GRAPH_DEBUG = 0
SAVE = 0

NORMALIZE = 0 # normalize 1 w.r.t. 2 (i.e. true w.r.t. reco)
              # 0: do not normalize
              # 1: normalize w.r.t. cube integral
              # 2: normalize w.r.t images integral (normalize each image on its own)

INPUT_DIR1 = '/home/mapaz/astropy/working_dir/bg_cube_models_comparison/gammapy_vs_hess-host/gammapy_scripts/20150831_for_hess_confluence/bg_cube_models_gammapy_a_la_michi'
BINNING_FORMAT1 = 'default'
NAME1 = 'gammapy'

INPUT_DIR2 = '/home/mapaz/HESS/fits_data/pa_fits_prod02/pa/Model_Deconvoluted_Prod26/Mpp_Std/background'
BINNING_FORMAT2 = 'michi'
NAME2 = 'michi'

E_REF = Quantity(1., 'TeV') # reference energy
NORM = 1
#INDEX = 2.7
#INDEX = 2.0
INDEX = 1.5

#SHOW_MODELS = True
SHOW_MODELS = False


def power_law(energy, E_0, norm, index):
    # TODO: use gammapy/spectrum/powerlaw.py for power-law functions!!!
    return norm*(energy/E_0)**-index


def int_power_law(energy_band, E_0, norm, index):
    # TODO: use gammapy/spectrum/powerlaw.py for power-law functions!!!
    return norm/(1 - index)/E_0**-index*(energy_band[1]**(1 - index) - energy_band[0]**(1 - index))


# group IDs for comparison

#  The following group IDs
#  - group_ids_selection = [14, 15, 20, 21, 26, 27]
# correspond to the cartesian product of the following
# "michi" alt az bin IDs:
#  - alt_bin_ids_selection = [7, 10, 13]
#  - az_bin_ids_selection = [0, 1]
group_ids_selection = [14, 15, 20, 21, 26, 27]


def convert_obs_groups_binning_def_michi_to_default():
    """Convert observation groups binning definition "michi" to "default".
    """
    # observation groups binning definition "michi"

    # alt az bin edges definitions
    altitude_edges = Angle([0, 20, 23, 27, 30, 33, 37, 40, 44, 49, 53, 58, 64, 72, 90], 'degree')
    azimuth_edges = Angle([-90, 90, 270], 'degree')

    # convert observation groups binning definition "michi" to "default"

    list_obs_group_axis = [ObservationGroupAxis('ALT', altitude_edges, 'bin_edges'),
                           ObservationGroupAxis('AZ', azimuth_edges, 'bin_edges')]
    obs_groups_michi = ObservationGroups(list_obs_group_axis)
    print("Observation groups 'michi':")
    print(obs_groups_michi.obs_groups_table)
    # save
    outfile = 'bg_observation_groups_michi.ecsv'
    print('Writing {}'.format(outfile))
    obs_groups_michi.write(outfile)

    # lookup table: equivalences in group/file naming "defualt" <-> "michi"
    # 3 columns: GROUP_ID, ALT_ID, AZ_ID
    # 28 rows: 1 per GROUP_ID

    lookup_obs_groups_michi = Table()
    n_cols = 1 + len(list_obs_group_axis)
    n_rows = obs_groups_michi.n_groups
    lookup_obs_groups_michi['GROUP_ID'] = np.zeros(n_rows, dtype=np.int)
    lookup_obs_groups_michi['ALT_ID'] = np.zeros(n_rows, dtype=np.int)
    lookup_obs_groups_michi['AZ_ID'] = np.zeros(n_rows, dtype=np.int)

    # loop over each observation group axis
    count_groups = 0
    for alt_id in np.arange(len(altitude_edges) - 1):
        for az_id in np.arange(len(azimuth_edges) - 1):
            lookup_obs_groups_michi['GROUP_ID'][count_groups] = count_groups
            lookup_obs_groups_michi['ALT_ID'][count_groups] = alt_id
            lookup_obs_groups_michi['AZ_ID'][count_groups] = az_id
            count_groups += 1

    print("lookup table:")
    print(lookup_obs_groups_michi)

    # save
    outfile = 'lookup_obs_groups_michi.ecsv'
    print('Writing {}'.format(outfile))
    # `~astropy.io.ascii` always overwrites the file
    ascii.write(lookup_obs_groups_michi, outfile,
                format='ecsv', fast_writer=False)


def look_obs_groups_michi(group_id):
    """
    Find corresponding ALT_ID, AZ_ID for a given GROUP_ID in lookup table.

    Parameters
    ----------
    group_id : int
        Group ID to look for.

    Returns
    -------
    i_alt : int
        Altitude bin index corresponding to the requested group ID.
    i_az : int
        Azimuth bin index corresponding to the requested group ID.
    """
    #read lookup
    filename = 'lookup_obs_groups_michi.ecsv'
    lookup_obs_groups_michi = ascii.read(filename)

    # find group row in lookup table
    group_ids = lookup_obs_groups_michi['GROUP_ID'].data
    group_index = np.where(group_ids == group_id)
    row = group_index[0][0]
    i_alt = lookup_obs_groups_michi['ALT_ID'][row]
    i_az = lookup_obs_groups_michi['AZ_ID'][row]
    return i_alt, i_az

# TODO: the lookup table is a bit hard-coded, but so is the binning
#       "michi"; otherwise I need to make a generic cartesian
#       product like the one I did in ObservationGroups for defining
#       the groups!!!
#       The same applies to the look function!!!


def plot_bg_cube_model_comparison(input_dir1, binning_format1, name1,
                                  input_dir2, binning_format2, name2):
    """
    Plot background cube model comparison.

    Produce a few figures for comparing 2 sets of bg cube models (1
    and 2), supposing same binning in both sets of observation
    groups (a.k.a. observation bin).

    Each figure corresponds to 1 observation group.
    Plot strategy in each figure:

    * Images:
        * rows: similar energy bin
        * cols: same bg cube model set
    * Spectra:
        * rows: similar det bin
        * cols: compare both bg cube model sets

    TODO: options to do more plots:

    * plot all groups:
        * leave 'group_ids_selection' empty
    * plot all bins:
        * USE 'CubeUtils' IN gammapy_bg_cube_models_utilities.py
          (TODO!!!)

    The script can be customized by setting a few global variables:

    * **group_ids_selection**: groups to compare; if empty: use all
      groups

    * **NORMALIZE**: normalization to use for the models. If
      activate, model 1 is normalized to match model 2. This can be
      useful, when comparing reco models w.r.t. true ones. Options:
          * *0*: do not normalize
          * *1*: normalize w.r.t. cube integral
          * *2*: normalize w.r.t images integral; each image (i.e.
            energy bin/slice) is normalized independently; this
            option can alter the spectral shape of the bg rate, but
            is is the way how smoothing method *michi* normalizes the
            background cube model, hence it is necessary to compare
            to those models that use that particular smoothing

    * **SAVE**: set to 1 (True) to save the output:
          * comparison plots as png
          * *michi* binning groups and lookup as ECVS files

    * **GRAPH_DEBUG**: if set to 1 (True) the program waits between
      each observation group iteration until the image is closed

    * **SHOW_MODELS**, **E_REF**, **NORM**, **INDEX**: to plot some
      power-laws normalized to the bg rate spectra for comparison

    Parameters
    ----------
    input_dir1, input_dir2 : str
        Directory where the corresponding set of bg cube models is stored.
    binning_format1, binning_format2 : {'default', 'michi'}
        String specifying the binning format; accepted values are:

        * *default* for the Gammapy format from
          `~gammapy.obs.ObservationGroups`; an observation groups
          ECVS file is expected in the bg cube models dir.
        * *michi* for the binning used by Michale Mayer;
          this script has methods to convert it to the
          *default* format.
          ref: [Mayer2015]_ (section 5.2.4)

    name1, name2 : str
        Name to use for plot labels/legends.

    References
    ==========
    .. [Mayer2015] `Michael Mayer (2015) <https://publishup.uni-potsdam.de/frontdoor/index/index/docId/7150>`_
       "Pulsar wind nebulae at high energies"
    """
    # check binning
    accepted_binnings = ['default', 'michi']

    if ((binning_format1 not in accepted_binnings) or
        (binning_format2 not in accepted_binnings)):
        raise ValueError("Invalid binning format: {0} or {1}".format(binning_format1,
                                                                     binning_format2))

    # convert binning, if necessary
    if binning_format1 == 'michi' or binning_format2 == 'michi':
        convert_obs_groups_binning_def_michi_to_default()

    # loop over observation groups: use binning of the 1st set to compare
    if binning_format1 == 'michi':
        observation_groups = obs_groups_michi
    else:
        observation_groups = ObservationGroups.read(input_dir1 + '/bg_observation_groups.ecsv')
    groups = observation_groups.list_of_groups
    print()
    print("list of groups", groups)

    for group in groups:
        print()
        print("group ", group)
        # compare only observation groups in group IDs selection
        # if empty, use all groups:
        if len(group_ids_selection) is not 0:
            groups_to_compare = group_ids_selection
        else:
            groups_to_compare = groups
        if group in groups_to_compare:
            group_info = observation_groups.info_group(group)
            print(group_info)

            # get cubes
            if binning_format1 == 'michi':
                # find corresponding ALT_ID, AZ_ID in lookup table
                i_alt, i_az = look_obs_groups_michi(group)
                filename1 = input_dir1 + '/hist_alt' + str(i_alt) +\
                            '_az' + str(i_az) + '.fits.gz'
            else:
                filename1 = input_dir1 + '/bg_cube_model_group' + str(group) +\
                            '_table.fits.gz'
            if binning_format2 == 'michi':
                # find corresponding ALT_ID, AZ_ID in lookup table
                i_alt, i_az = look_obs_groups_michi(group)
                filename2 = input_dir2 + '/hist_alt' + str(i_alt) +\
                            '_az' + str(i_az) + '.fits.gz'
            else:
                filename2 = input_dir2 + '/bg_cube_model_group' + str(group) +\
                            '_table.fits.gz'
            print('filename1', filename1)
            print('filename2', filename2)
            bg_cube_model1 = CubeBackgroundModel.read(filename1,
                                                      format='table').background_cube
            bg_cube_model2 = CubeBackgroundModel.read(filename2,
                                                      format='table').background_cube
            #bg_cube_model2 = CubeBackgroundModel.read(filename2,
            #                                          format='table').counts_cube # events!!!

            # normalize 1 w.r.t. 2 (i.e. true w.r.t. reco)
            if NORMALIZE == 1:
                # normalize w.r.t. cube integral
                integral1 = bg_cube_model1.integral
                integral2 = bg_cube_model2.integral
                bg_cube_model1.data *= integral2/integral1
            elif NORMALIZE == 2:
                # normalize w.r.t images integral (normalize each image on its own)
                integral_images1 = bg_cube_model1.integral_images
                integral_images2 = bg_cube_model2.integral_images
                for i_energy in np.arange(len(bg_cube_model1.energy_edges) - 1):
                    bg_cube_model1.data[i_energy] *= (integral_images2/integral_images1)[i_energy]

            # compare binning
            print("energy edges 1", bg_cube_model1.energy_edges)
            print("energy edges 2", bg_cube_model2.energy_edges)
            print("detector edges 1 Y", bg_cube_model1.coordy_edges)
            print("detector edges 2 Y", bg_cube_model2.coordy_edges)
            print("detector edges 1 X", bg_cube_model1.coordx_edges)
            print("detector edges 2 X", bg_cube_model2.coordx_edges)

            # make sure that both cubes use the same units for the plots
            bg_cube_model2.data = bg_cube_model2.data.to(bg_cube_model1.data.unit)
            #bg_cube_model1.data = bg_cube_model1.data.to('1 / (TeV s sr)')
            #bg_cube_model2.data = bg_cube_model2.data.to('1 / (TeV s sr)')
            # TODO: actually the display is better if using '1 / (TeV s sr)'

            # plot
            fig, axes = plt.subplots(nrows=2, ncols=3)
            fig.set_size_inches(30., 15., forward=True)
            plt.suptitle(group_info)

            # plot images
            #  rows: similar energy bin
            #  cols: same file
            #bg_cube_model1.plot_image(energy=Quantity(0.5, 'TeV'), ax=axes[0, 0])
            bg_cube_model1.plot_image(energy=Quantity(5., 'TeV'), ax=axes[0, 0])
            axes[0, 0].set_title("{0}: {1}".format(name1, axes[0, 0].get_title()))
            bg_cube_model1.plot_image(energy=Quantity(50., 'TeV'), ax=axes[1, 0])
            axes[1, 0].set_title("{0}: {1}".format(name1, axes[1, 0].get_title()))
            #bg_cube_model2.plot_image(energy=Quantity(0.5, 'TeV'), ax=axes[0, 1])
            bg_cube_model2.plot_image(energy=Quantity(5., 'TeV'), ax=axes[0, 1])
            axes[0, 1].set_title("{0}: {1}".format(name2, axes[0, 1].get_title()))
            bg_cube_model2.plot_image(energy=Quantity(50., 'TeV'), ax=axes[1, 1])
            axes[1, 1].set_title("{0}: {1}".format(name2, axes[1, 1].get_title()))

            # plot spectra
            #  rows: similar det bin
            #  cols: compare both files
            bg_cube_model1.plot_spectrum(coord=Angle([0., 0.], 'degree'),
                                         ax=axes[0, 2],
                                         style_kwargs=dict(color='blue',
                                                           label=name1))
            spec_title1 = axes[0, 2].get_title()
            bg_cube_model2.plot_spectrum(coord=Angle([0., 0.], 'degree'),
                                         ax=axes[0, 2],
                                         style_kwargs=dict(color='red',
                                                           label=name2))
            spec_title2 = axes[0, 2].get_title()
            if spec_title1 != spec_title2:
                s_error = "Expected same det binning, but got "
                s_error += "\"{0}\" and \"{1}\"".format(spec_title1, spec_title2)
                raise ValueError(s_error)
            else:
                axes[0, 2].set_title(spec_title1)

            if SHOW_MODELS:
                # plot normalized models on top

                E_0 = E_REF
                norm = NORM
                index = INDEX

                plot_data_x = axes[0, 2].get_lines()[0].get_xydata()[:,0]
                plot_data_y = axes[0, 2].get_lines()[0].get_xydata()[:,1]
                plot_data_int = np.trapz(y=plot_data_y, x=plot_data_x)
                energy_band = np.array([plot_data_x[0], plot_data_x[-1]])
                model_int = int_power_law(energy_band, E_0, norm, index)
                normed_PL1 = plot_data_int/model_int*power_law(plot_data_x, E_0, norm, index)
                axes[0, 2].plot(plot_data_x, normed_PL1, color='blue',
                                linestyle='dotted', linewidth=2,
                                label='model index = {}'.format(index))

                index = INDEX + 1

                plot_data_x = axes[0, 2].get_lines()[0].get_xydata()[:,0]
                plot_data_y = axes[0, 2].get_lines()[0].get_xydata()[:,1]
                plot_data_int = np.trapz(y=plot_data_y, x=plot_data_x)
                energy_band = np.array([plot_data_x[0], plot_data_x[-1]])
                model_int = int_power_law(energy_band, E_0, norm, index)
                normed_PL2 = plot_data_int/model_int*power_law(plot_data_x, E_0, norm, index)
                axes[0, 2].plot(plot_data_x, normed_PL2, color='blue',
                                linestyle='dashed', linewidth=2,
                                label='model index = {}'.format(index))

            axes[0, 2].legend()

            bg_cube_model1.plot_spectrum(coord=Angle([2., 2.], 'degree'),
                                         ax=axes[1, 2],
                                         style_kwargs=dict(color='blue',
                                                           label=name1))
            spec_title1 = axes[1, 2].get_title()
            bg_cube_model2.plot_spectrum(coord=Angle([2., 2.], 'degree'),
                                         ax=axes[1, 2],
                                         style_kwargs=dict(color='red',
                                                           label=name2))
            spec_title2 = axes[1, 2].get_title()
            if spec_title1 != spec_title2:
                s_error = "Expected same det binning, but got "
                s_error += "\"{0}\" and \"{1}\"".format(spec_title1, spec_title2)
                raise ValueError(s_error)
            else:
                axes[1, 2].set_title(spec_title1)

            if SHOW_MODELS:
                # plot normalized models on top

                E_0 = E_REF
                norm = NORM
                index = INDEX

                plot_data_x = axes[1, 2].get_lines()[0].get_xydata()[:,0]
                plot_data_y = axes[1, 2].get_lines()[0].get_xydata()[:,1]
                plot_data_int = np.trapz(y=plot_data_y, x=plot_data_x)
                energy_band = np.array([plot_data_x[0], plot_data_x[-1]])
                model_int = int_power_law(energy_band, E_0, norm, index)
                normed_PL1 = plot_data_int/model_int*power_law(plot_data_x, E_0, norm, index)
                axes[1, 2].plot(plot_data_x, normed_PL1, color='blue',
                                linestyle='dotted', linewidth=2,
                                label='model index = {}'.format(index))

                index = INDEX + 1

                plot_data_x = axes[1, 2].get_lines()[0].get_xydata()[:,0]
                plot_data_y = axes[1, 2].get_lines()[0].get_xydata()[:,1]
                plot_data_int = np.trapz(y=plot_data_y, x=plot_data_x)
                energy_band = np.array([plot_data_x[0], plot_data_x[-1]])
                model_int = int_power_law(energy_band, E_0, norm, index)
                normed_PL2 = plot_data_int/model_int*power_law(plot_data_x, E_0, norm, index)
                axes[1, 2].plot(plot_data_x, normed_PL2, color='blue',
                                linestyle='dashed', linewidth=2,
                                label='model index = {}'.format(index))

            axes[1, 2].legend()

            if GRAPH_DEBUG:
                plt.show() # wait until image is closed

            if SAVE:
                outfile = "bg_cube_model_comparison_group{}.png".format(group)
                print('Writing {}'.format(outfile))
                fig.savefig(outfile)

    plt.show() # don't leave at the end


if __name__ == '__main__':
    """Main function: parse arguments and launch the whole analysis chain.
    """
    parser = argparse.ArgumentParser(description='Compare 2 sets of bg cube models.')

    parser.add_argument('--input-dir1', type=str,
                        default=INPUT_DIR1,
                        help='Directory where the corresponding set '
                        'of bg cube models is stored. '
                        '(default is {}).'.format(INPUT_DIR1))

    parser.add_argument('--binning-format1', type=str,
                        default=BINNING_FORMAT1,
                        choices=['default', 'michi'],
                        help='String specifying the binning format. '
                        '(default is {}).'.format(BINNING_FORMAT1))

    parser.add_argument('--name1', type=str,
                        default=NAME1,
                        help='Name to use for plot labels/legends. '
                        '(default is {}).'.format(NAME1))

    parser.add_argument('--input-dir2', type=str,
                        default=INPUT_DIR2,
                        help='Directory where the corresponding set '
                        'of bg cube models is stored. '
                        '(default is {}).'.format(INPUT_DIR2))

    parser.add_argument('--binning-format2', type=str,
                        default=BINNING_FORMAT2,
                        choices=['default', 'michi'],
                        help='String specifying the binning format. '
                        '(default is {}).'.format(BINNING_FORMAT2))

    parser.add_argument('--name2', type=str,
                        default=NAME2,
                        help='Name to use for plot labels/legends. '
                        '(default is {}).'.format(NAME2))

    args = parser.parse_args()

    plot_bg_cube_model_comparison(args.input_dir1, args.binning_format1, args.name1,
                                  args.input_dir2, args.binning_format2, args.name2)

from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
import numpy as np
import matplotlib.pyplot as plt
from astropy.units import Quantity, UnitsError
from astropy.coordinates import Angle
from astropy.table import Table
from gammapy.obs import DataStore
from gammapy.data import EventListDataset

GRAPH_DEBUG = 0
SAVE = 0

FITSPATH = '/home/mapaz/astropy/working_dir/bg_cube_models_comparison/gammapy_true_vs_reco/20150901_debugging_true_reco' + '/test_dataset'
SCHEME = 'HESS'

# model parameters: TODO: match them to the dataset!!!

E_REF = Quantity(1., 'TeV') # reference energy
#NORM = Quantity(1., '1 / (s TeV sr)')
NORM = 1
#INDEX = 2.7
#INDEX = 2.0
INDEX = 1.5


def power_law(energy, E_0, norm, index):
    # TODO: use gammapy/spectrum/powerlaw.py for power-law functions!!!
    return norm*(energy/E_0)**-index


def int_power_law(energy_band, E_0, norm, index):
    # TODO: use gammapy/spectrum/powerlaw.py for power-law functions!!!
    return norm/(1 - index)/E_0**-index*(energy_band[1]**(1 - index) - energy_band[0]**(1 - index))


def plot_dataset():
    fits_path = FITSPATH
    scheme = SCHEME
    data_store = DataStore(dir=fits_path, scheme=scheme)
    observation_table = data_store.make_observation_table()
    event_list_files = data_store.make_table_of_files(observation_table,
                                                          filetypes=['events'])
    data_set = EventListDataset.vstack_from_files(event_list_files['filename']).event_list
    # the stacking can be long, if many runs are read: if only a few
    # variables are needed it is faster to grab the needed columns
    # and stack them manually (like for the
    # CubeBackgroundModel.fill_events method)

    # construct events table (X, Y, energy)
    # TODO: units are missing in the H.E.S.S. fits event
    #       lists; this should be solved in the next (prod03)
    #       H.E.S.S. fits production
    # workaround: try to cast units, if it doesn't work, use hard coded
    # ones
    try:
        ev_DETX = Angle(data_set['DETX'])
        ev_DETY = Angle(data_set['DETY'])
        ev_energy = Quantity(data_set['ENERGY'])
    except UnitsError:
        ev_DETX = Angle(data_set['DETX'], 'degree') # hard-coded!!!
        ev_DETY = Angle(data_set['DETY'], 'degree') # hard-coded!!!
        ev_energy = Quantity(data_set['ENERGY'],
                             data_set.meta['EUNIT']) # half hard-coded!!!
    ev_table = Table([ev_DETX, ev_DETY, ev_energy],
                          names=('DETX', 'DETY', 'ENERGY'))
    print(ev_table)

    ############
    # plotting #
    ############

    # model parameters

    E_0 = E_REF
    norm = NORM
    index = INDEX
    #index = index - 1 ## seems to work for sample_PL (it does exp = 1 - gamma, why??!!!) -> workaround: use sample_PL(index+1)

    n_events = len(ev_table) # not used in the current version!!!

    # energy histogram: log binning

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # energy bins (logarithmic)
    #energy_band = Quantity([0.1, 80.], 'TeV')
    energy_band = Quantity([0.1, 100.], 'TeV')
    nenergy_bins = 20
    log_delta_energy = (np.log(energy_band[1].value)
                        - np.log(energy_band[0].value))/nenergy_bins
    energy_edges = np.exp(np.arange(nenergy_bins + 1)*log_delta_energy
                          + np.log(energy_band[0].value))
    energy_edges = Quantity(energy_edges, energy_band[0].unit)
    count, bins, ignored = ax.hist(ev_energy.value, energy_edges.value)

    # plot normalized model on top
    hist_int = (count*np.diff(bins)).sum()
    model_int = int_power_law(energy_band, E_0, norm, index)
    normed_PL = hist_int/model_int*power_law(bins, E_0, norm, index)
    ax.plot(bins, normed_PL, color='red')

    ax.loglog() # double log scale # slow!

    if GRAPH_DEBUG:
        plt.show() # wait until image is closed

    if SAVE:
        outfile = "dataset_energy_hist.png"
        print('Writing {}'.format(outfile))
        fig.savefig(outfile)

    plt.show() #don't quit at the end


if __name__ == '__main__':
    plot_dataset()

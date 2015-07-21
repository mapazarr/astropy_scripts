from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
from astropy.units import Quantity
from astropy.coordinates import Angle
from gammapy.background import CubeBackgroundModel

__all__ = ['CubeBackgroundModelUtils',
           ]

class CubeBackgroundModelUtils():

    def plot_images(bg_model, energy=None, style_kwargs=None):
        """Plot images for each energy bin.

        Save images in files: several canvases with a few images
        each, 1 file per canvas.
        If specifying a particular energy, the function returns the
        figure of the specific energy bin containing the specified
        value. If no energy is specified, no figure is returned,
        since it would be very memory consuming.

        TODO: solve the problem of having to choose between: (!!!)

        * Linear color axis: showing horrible non scientific notation
          ticks (i.e. 0.0001).
        * Log color axis (LogNorm()): showing all images almost equal
          to each other.

        ref: http://stackoverflow.com/questions/25983218/scientific-notation-colorbar-in-matplotlib (!!!)

        Parameters
        ----------
        bg_model : `~gammapy.background.CubeBackgroundModel`
            cube bg model object to plot
        energy : `~astropy.units.Quantity`, optional
                energy of bin to plot the bg model
        style_kwargs : `~dict`
                style arguments to pass to the plotting function

        Returns
        -------
        axes : `~matplotlib.axes.Axes`
            axes of the figure with image of bin of the bg model
            for the selected energy value (if any), optional
        """
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        do_only_1_plot = False # general case: print all plots
        if energy:
            energy = energy.flatten() # flatten
            # check shape of energy: only 1 value is accepted
            nvalues = len(energy)
            if nvalues != 1:
                raise IndexError("Expected exactly 1 value for energy, got {}."
                                 .format(nvalues))
            else:
                energy = Quantity(energy[0])
                print("Reqested plot only for 1 energy: {}".format(energy))
                do_only_1_plot = True

        n_energy_bins = len(bg_model.energy_bins) - 1
        nimages = n_energy_bins
        ncols = 4
        nrows = 4
        if do_only_1_plot:
            n_energy_bins = nimages = ncols = nrows = 1
        npads_per_canvas = ncols*nrows

        fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
        fig.set_size_inches(35., 35., forward=True)
        count_images = 1
        count_canvases = 1
        count_pads = 1

        if style_kwargs == None:
            style_kwargs=dict()

        extent = bg_model.image_extent
        energy_bin_centers = bg_model.energy_bin_centers
        if do_only_1_plot:
            # find energy bin containing the specified energy
            energy_bin = bg_model.find_energy_bin(energy)
            energy_bin_edges = bg_model.find_energy_bin_edges(energy)
            ss_energy_bin_edges = "[{0}, {1}) {2}".format(energy_bin_edges[0].value,
                                                          energy_bin_edges[1].value,
                                                          energy_bin_edges.unit)
            print("Found energy {0} in bin {1} with boundaries {2}.".format(energy,
                                                                            energy_bin,
                                                                            ss_energy_bin_edges))

        for ii in range(n_energy_bins):
            if do_only_1_plot:
                ii = energy_bin
            print("ii", ii)
            data = bg_model.background[ii]
            energy_bin_center = energy_bin_centers[ii]
            print ("  image({}) canvas({}) pad({})".format(count_images, count_canvases, count_pads))

            if do_only_1_plot:
                fig.set_size_inches(8., 8., forward=True)
                ax = axes
            else:
                ax = axes.flat[count_pads - 1]
            ax = bg_model.plot_image(energy_bin_center, ax, style_kwargs)

            if do_only_1_plot:
                ax.set_title('Energy = [{0:.1f}, {1:.1f}) {2}'.format(energy_bin_edges[0].value,
                                                                      energy_bin_edges[1].value,
                                                                      energy_bin_edges.unit))
            else:
                ax.set_title('Energy = {:.1f}'.format(energy_bin_center))

            count_pads += 1 # increase

            if count_pads > npads_per_canvas or count_images == nimages:
                print("Canvas full, saving and creating a new canvas")
                if do_only_1_plot:
                    filename = "cube_background_model_image{0}{1}.png".format(energy.value, energy.unit)
                else:
                    filename = "cube_background_model_images{}.png".format(count_canvases)
                print('Writing {}'.format(filename))
                fig.savefig(filename)
                if not do_only_1_plot:
                    plt.close('all') # close all open figures
                    fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
                    fig.set_size_inches(35., 35., forward=True)
                count_canvases += 1 # increase
                count_pads = 1 # reset

            count_images += 1 # increase

        if do_only_1_plot:
            return axes

    def plot_spectra(bg_model, format=None, det=None, style_kwargs=None):
        """Plot spectra for each spatial (X, Y) bin.

        Save images in files: several canvases with a few images
        each, 1 file per canvas.
        If specifying a particular det (X,Y) pair, the function
        returns the figure of the specific det bin containing the
        specified value. If no det is specified, no figure is
        returned, since it would be very memory consuming.

        Several plotting formats are accepted, depending on the value
        of the `format` parameter:

            * `mosaic`: mosaic of figures with several pads per figure

            * `stack`: 1 plot with several curves on the same axis

        Parameters
        ----------
        bg_model : `~gammapy.background.CubeBackgroundModel`
            cube bg model object to plot
        format : `~string`, optional
            format for the plots
        det : `~astropy.units.Quantity`, optional
            det (X,Y) pair of bin to plot the bg model
        style_kwargs : `~dict`
                style arguments to pass to the plotting function

        Returns
        -------
        axes : `~matplotlib.axes.Axes`
            axes of the figure with image of bin of the bg model
            for the selected det (X,Y) pair (if any), optional
        """
        import matplotlib.pyplot as plt

        do_only_1_plot = False # general case: print all plots
        if det:
            det = det.flatten() # flatten
            # check shape of det: only 1 pair is accepted
            nvalues = len(det.flatten())
            if nvalues != 2:
                ss_error = "Expected exactly 2 values for det (X, Y),"
                ss_error += "got {}.".format(nvalues)
                raise IndexError(ss_error)
            else:
                print("Reqested plot only for 1 det: {}".format(det))
                do_only_1_plot = True

        # check format
        valid_formats = ['mosaic', 'stack']
        if do_only_1_plot:
            format = None
        else:
            if format not in valid_formats:
                raise ValueError("Plot format {} not understood!".format(format))

        n_det_bins_x = len(bg_model.detx_bins) - 1
        n_det_bins_y = len(bg_model.dety_bins) - 1
        nimages = n_det_bins_x*n_det_bins_y
        ncols = 4
        nrows = 4
        if do_only_1_plot:
            n_det_bins = n_det_bins_x = n_det_bins_y = nimages = ncols = nrows = 1
        npads_per_canvas = ncols*nrows
        if format == 'stack':
            # reset ncols and nrows
            # npads_per_canvas will be used as number of curves per
            # canvas (or pad, since there is only 1 pad per canvas)
            ncols = nrows = 1

        fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
        fig.set_size_inches(25., 25., forward=True)
        count_images = 1
        count_canvases = 1
        count_pads = 1

        if style_kwargs == None:
            style_kwargs=dict()

        energy_points = bg_model.energy_bin_centers
        detx_bin_centers, dety_bin_centers = bg_model.image_bin_centers
        if do_only_1_plot:
            # find det bin containing the specified det coordinates
            det_bin = bg_model.find_det_bin(det)
            det_bin_edges = bg_model.find_det_bin_edges(det)
            ss_detx_bin_edges = "[{0}, {1}) {2}".format(det_bin_edges[0].value,
                                                        det_bin_edges[1].value,
                                                        det_bin_edges.unit)
            ss_dety_bin_edges = "[{0}, {1}) {2}".format(det_bin_edges[2].value,
                                                        det_bin_edges[3].value,
                                                        det_bin_edges.unit)
            print("Found det {0} in bin {1} with boundaries {2}, {3}.".format(det, det_bin,
                                                                              ss_detx_bin_edges,
                                                                              ss_dety_bin_edges))

        for ii in range(n_det_bins_x):
            if do_only_1_plot:
                ii = det_bin[0]
            print("ii", ii)
            for jj in range(n_det_bins_y):
                if do_only_1_plot:
                    jj = det_bin[1]
                print(" jj", jj)
                data = bg_model.background[:, jj, ii]
                detx_bin_center = detx_bin_centers[ii]
                dety_bin_center = dety_bin_centers[jj]
                det_bin_center = Angle([detx_bin_center, dety_bin_center])
                if format == 'stack':
                    print ("  image({}) canvas({}) plot({})".format(count_images, count_canvases, count_pads))
                else:
                    print ("  image({}) canvas({}) pad({})".format(count_images, count_canvases, count_pads))

                if do_only_1_plot or format == 'stack':
                    fig.set_size_inches(8., 8., forward=True)
                    ax = axes
                else:
                    ax = axes.flat[count_pads - 1]
                ss_det_bin_center = "({0:.1f}, {1:.1f})".format(detx_bin_center, dety_bin_center)
                ss_label = 'Det = {}'.format(ss_det_bin_center)
                style_kwargs['label'] = ss_label
                ax = bg_model.plot_spectrum(det_bin_center, ax, style_kwargs)
                if do_only_1_plot:
                    ss_detx_bin_edges = "[{0:.1f}, {1:.1f}) {2}".format(det_bin_edges[0].value,
                                                                        det_bin_edges[1].value,
                                                                        det_bin_edges.unit)
                    ss_dety_bin_edges = "[{0:.1f}, {1:.1f}) {2}".format(det_bin_edges[2].value,
                                                                        det_bin_edges[3].value,
                                                                        det_bin_edges.unit)

                    ax.set_title('Det = {0} {1}'.format(ss_detx_bin_edges, ss_dety_bin_edges))
                else:
                    if format == 'stack':
                        ax.set_title('')
                    else:
                        ax.set_title(ss_label)
                count_pads += 1 # increase

                if count_pads > npads_per_canvas or count_images == nimages:
                    print("Canvas full, saving and creating a new canvas")
                    if do_only_1_plot:
                        filename = "cube_background_model_spectrum{0}{2}{1}{2}.png".format(det.value[0],
                                                                                           det.value[1],
                                                                                           det.unit)
                    else:
                        if format == 'stack':
                            ax.legend()
                        filename = "cube_background_model_spectra_{0}{1}.png".format(format, count_canvases)
                    print('Writing {}'.format(filename))
                    fig.savefig(filename)
                    if not do_only_1_plot:
                        plt.close('all') # close all open figures
                        fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
                        if format == 'stack':
                            fig.set_size_inches(8., 8., forward=True)
                        else:
                            fig.set_size_inches(25., 25., forward=True)
                    count_canvases += 1 # increase
                    count_pads = 1 # reset

                count_images += 1 # increase

        if do_only_1_plot:
            return axes

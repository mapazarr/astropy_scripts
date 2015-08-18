from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
from astropy.units import Quantity
from astropy.coordinates import Angle
from gammapy.background import Cube

__all__ = ['CubeUtils',
           ]

class CubeUtils():

    def plot_images(cube, energy=None, style_kwargs=None):
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

        TODO: generalize naming scheme for all kinds of cubes!!!
        TODO: rename this file to gammapy_cube_utilities.py !!! (and update references to it in other scripts)!!!

        Parameters
        ----------
        cube : `~gammapy.background.Cube`
            cube object to plot
        energy : `~astropy.units.Quantity`, optional
                energy of bin to plot the cube
        style_kwargs : `~dict`
                style arguments to pass to the plotting function

        Returns
        -------
        axes : `~matplotlib.axes.Axes`
            axes of the figure with image of bin of the cube
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

        n_energy_bins = len(cube.energy_edges) - 1
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

        extent = cube.image_extent
        energy_bin_centers = cube.energy_bin_centers
        if do_only_1_plot:
            # find energy bin containing the specified energy
            energy_bin = cube.find_energy_bin(energy)
            energy_bin_edges = cube.find_energy_bin_edges(energy)
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
            data = cube.data[ii]
            energy_bin_center = energy_bin_centers[ii]
            print ("  image({}) canvas({}) pad({})".format(count_images, count_canvases, count_pads))

            if do_only_1_plot:
                fig.set_size_inches(8., 8., forward=True)
                ax = axes
            else:
                ax = axes.flat[count_pads - 1]
            ax = cube.plot_image(energy_bin_center, ax, style_kwargs)

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

    def plot_spectra(cube, format=None, coord=None, style_kwargs=None):
        """Plot spectra for each spatial (X, Y) bin.

        Save images in files: several canvases with a few images
        each, 1 file per canvas.
        If specifying a particular coord (X,Y) pair, the function
        returns the figure of the specific coord bin containing the
        specified value. If no coord is specified, no figure is
        returned, since it would be very memory consuming.

        Several plotting formats are accepted, depending on the value
        of the `format` parameter:

            * `mosaic`: mosaic of figures with several pads per figure

            * `stack`: 1 plot with several curves on the same axis

        Parameters
        ----------
        cube : `~gammapy.background.Cube`
            cube object to plot
        format : `~string`, optional
            format for the plots
        coord : `~astropy.units.Quantity`, optional
            coord (X,Y) pair of bin to plot the cube
        style_kwargs : `~dict`
                style arguments to pass to the plotting function

        Returns
        -------
        axes : `~matplotlib.axes.Axes`
            axes of the figure with image of bin of the cube
            for the selected coord (X,Y) pair (if any), optional
        """
        import matplotlib.pyplot as plt

        do_only_1_plot = False # general case: print all plots
        if coord:
            coord = coord.flatten() # flatten
            # check shape of coord: only 1 pair is accepted
            nvalues = len(coord.flatten())
            if nvalues != 2:
                ss_error = "Expected exactly 2 values for coord (X, Y),"
                ss_error += "got {}.".format(nvalues)
                raise IndexError(ss_error)
            else:
                print("Reqested plot only for 1 coord: {}".format(coord))
                do_only_1_plot = True

        # check format
        valid_formats = ['mosaic', 'stack']
        if do_only_1_plot:
            format = None
        else:
            if format not in valid_formats:
                raise ValueError("Plot format {} not understood!".format(format))

        n_coord_bins_x = len(cube.coordx_edges) - 1
        n_coord_bins_y = len(cube.coordy_edges) - 1
        nimages = n_coord_bins_x*n_coord_bins_y
        ncols = 4
        nrows = 4
        if do_only_1_plot:
            n_coord_bins = n_coord_bins_x = n_coord_bins_y = nimages = ncols = nrows = 1
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

        energy_points = cube.energy_bin_centers
        coordx_bin_centers, coordy_bin_centers = cube.image_bin_centers
        if do_only_1_plot:
            # find coord bin containing the specified coord coordinates
            coord_bin = cube.find_coord_bin(coord)
            coord_bin_edges = cube.find_coord_bin_edges(coord)
            ss_coordx_bin_edges = "[{0}, {1}) {2}".format(coord_bin_edges[0].value,
                                                        coord_bin_edges[1].value,
                                                        coord_bin_edges.unit)
            ss_coordy_bin_edges = "[{0}, {1}) {2}".format(coord_bin_edges[2].value,
                                                        coord_bin_edges[3].value,
                                                        coord_bin_edges.unit)
            print("Found coord {0} in bin {1} with boundaries {2}, {3}.".format(coord, coord_bin,
                                                                              ss_coordx_bin_edges,
                                                                              ss_coordy_bin_edges))

        for ii in range(n_coord_bins_x):
            if do_only_1_plot:
                ii = coord_bin[0]
            print("ii", ii)
            for jj in range(n_coord_bins_y):
                if do_only_1_plot:
                    jj = coord_bin[1]
                print(" jj", jj)
                data = cube.data[:, jj, ii]
                coordx_bin_center = coordx_bin_centers[ii]
                coordy_bin_center = coordy_bin_centers[jj]
                coord_bin_center = Angle([coordx_bin_center, coordy_bin_center])
                if format == 'stack':
                    print ("  image({}) canvas({}) plot({})".format(count_images, count_canvases, count_pads))
                else:
                    print ("  image({}) canvas({}) pad({})".format(count_images, count_canvases, count_pads))

                if do_only_1_plot or format == 'stack':
                    fig.set_size_inches(8., 8., forward=True)
                    ax = axes
                else:
                    ax = axes.flat[count_pads - 1]
                ss_coord_bin_center = "({0:.1f}, {1:.1f})".format(coordx_bin_center, coordy_bin_center)
                ss_label = 'Coord = {}'.format(ss_coord_bin_center)
                style_kwargs['label'] = ss_label
                ax = cube.plot_spectrum(coord_bin_center, ax, style_kwargs)
                if do_only_1_plot:
                    ss_coordx_bin_edges = "[{0:.1f}, {1:.1f}) {2}".format(coord_bin_edges[0].value,
                                                                        coord_bin_edges[1].value,
                                                                        coord_bin_edges.unit)
                    ss_coordy_bin_edges = "[{0:.1f}, {1:.1f}) {2}".format(coord_bin_edges[2].value,
                                                                        coord_bin_edges[3].value,
                                                                        coord_bin_edges.unit)

                    ax.set_title('Coord = {0} {1}'.format(ss_coordx_bin_edges, ss_coordy_bin_edges))
                else:
                    if format == 'stack':
                        ax.set_title('')
                    else:
                        ax.set_title(ss_label)
                count_pads += 1 # increase

                if count_pads > npads_per_canvas or count_images == nimages:
                    print("Canvas full, saving and creating a new canvas")
                    if do_only_1_plot:
                        filename = "cube_background_model_spectrum{0}{2}{1}{2}.png".format(coord.value[0],
                                                                                           coord.value[1],
                                                                                           coord.unit)
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

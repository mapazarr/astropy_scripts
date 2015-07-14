import numpy as np
import matplotlib.pyplot as plt
from plot_something import plot_something

GRAPH_DEBUG = 1

# example 1: 2 plots in the same axis

# x coord
x = np.arange(0, 10, 0.1)

# 1st plot
y = np.sin(x)
ax = plot_something(x, y, style_kwargs=dict(color='blue', label='sin'))

# 2nd plot on top (like Draw("same") in ROOT)
y2 = np.cos(x)
ax = plot_something(x, y2, ax=ax, style_kwargs=dict(color='red', label='cos'))

# legend
ax.legend()

plt.draw()
if GRAPH_DEBUG:
    plt.show() # wait until image is closed

# example 2: 2 plots in 2 pads of the same canvas

fig, axes = plt.subplots(nrows=1, ncols=2)
fig.set_size_inches(16., 8., forward=True)

axes[0] = plot_something(x, y, ax=axes[0], style_kwargs=dict(color='blue', label='sin'))
axes[0].legend()

axes[1] = plot_something(x, y2, ax=axes[1], style_kwargs=dict(color='red', label='cos'))
axes[1].legend()

plt.draw()
if GRAPH_DEBUG:
    plt.show() # wait until image is closed

# examples using CubeBackgroundModel class plots

from astropy.units import Quantity
from astropy.coordinates import Angle
from gammapy.background import CubeBackgroundModel
from gammapy import datasets

# example 3: same as example 2, but using plot_image function
#            from CubeBackgroundModel class

filename = '../test_datasets/background/bg_cube_model_test.fits'
filename = datasets.get_path(filename, location='remote')
bg_model = CubeBackgroundModel.read(filename, format='bin_table')

fig, axes = plt.subplots(nrows=1, ncols=2)
fig.set_size_inches(16., 8., forward=True)

axes[0] = bg_model.plot_image(energy=Quantity(2., 'TeV'), ax=axes[0])
axes[1] = bg_model.plot_image(energy=Quantity(20., 'TeV'), ax=axes[1])

plt.draw()
if GRAPH_DEBUG:
    plt.show() # wait until image is closed

# example 4: same as example 1, but using plot_spectrum function
#            from CubeBackgroundModel class





# example 5: plot all images (i.e. one image per energy slice) in bg cube





# example 6: plot all spectra (i.e. one spectrum per det (X, Y) bin) in bg cube






plt.show() #don't quit at the end

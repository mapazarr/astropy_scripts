from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm
from pylab import *
import matplotlib.mlab as mlab

GRAPH_DEBUG = 0

#################
# histogramming #
#################

data = np.array([1.,2.,3.,4.,5.])
print("data", repr(data))
bins = np.array([0., 2.5, 5.])
print("bins", repr(bins))
dig_data = np.digitize(data, bins)
h_data, edges_data = np.histogram(data, bins)

data_x = np.array([1.,2.,3.,4.,5.]) # column-wise # works for histogram2d
data_y = np.array([11.,12.,13.,14.,15.]) # works for histogram2d
print("data_x", repr(data_x))
print("data_y", repr(data_y))
data_2d = np.array([[1.,2.,3.,4.,5.], [11.,12.,13.,14.,15.]]) # column-wise # doesn't work!
data_2d_col = np.array([[1.,2.,3.,4.,5.], [11.,12.,13.,14.,15.]]) # column-wise # doesn't work!
data_2d_row = np.array([[1.,11.], [2.,12.], [3.,13.], [4.,14.], [5.,15.]]) # row-wise # works for histogramdd
print("data_2d", repr(data_2d))
#bins_2d = np.array([[0., 2.5, 5.], [10., 12.5, 15.]])
bins_2d = np.array([[0., 2., 4., 6.], [10., 12., 14., 16.]])
print("bins_2d", repr(bins_2d))
#np.digitize(data_2d, bins_2d) # doesn't work!
#np.histogram(data_2d, bins_2d) # doesn't work!
#np.histogram2d(data_2d, bins_2d) # doesn't work!
h2d_data, x_edges_data, y_edges_data = np.histogram2d(data_x, data_y, bins_2d)
hdd_data, edges_data = np.histogramdd(data_2d_row, bins_2d)

# convert column-wise data to row-wise data
data_2d_col_trans = np.transpose(data_2d_col)

# test histogram from astropy table
astro_table = Table()
astro_table['X'] = data_x
astro_table['Y'] = data_y
print("astro_table")
print(astro_table)
data = np.array(astro_table) # structured array incompatible with histogram functions!
print("structured data", data)
data = np.vstack([astro_table['X'], astro_table['Y']]).T
print("non-structured data")
print(data)
hdd_table, edges_table = np.histogramdd(data, bins_2d)

#import IPython; IPython.embed()
#dig_data
#h_data
#h2d_data
#hdd_data
#hdd_table
print("dig_data")
print(repr(dig_data))
print("h_data")
print(repr(h_data))
print("h2d_data")
print(repr(h2d_data))
print("hdd_data")
print(repr(hdd_data))
print("hdd_table")
print(repr(hdd_table))

# there is also methods to histogram directly when plotting with mpl (plt.hist)
# ref: http://matplotlib.org/examples/statistics/histogram_demo_features.html

print()

############
# plotting #
############

# Examples of 1D histogram plots
# ref: http://matplotlib.org/examples/statistics/histogram_demo_features.html
#      http://matplotlib.org/examples/statistics/histogram_demo_histtypes.html

# 1st example

fig = plt.figure()
ax = fig.add_subplot(111)

# example data
mu = 100 # mean of distribution
sigma = 15 # standard deviation of distribution
x = mu + sigma * np.random.randn(10000)

num_bins = 50
# the histogram of the data
n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
# add a 'best fit' line
y = mlab.normpdf(bins, mu, sigma)
plt.plot(bins, y, 'r--')
plt.xlabel('Smarts')
plt.ylabel('Probability')
plt.title(r'Histogram of IQ: $\mu=100$, $\sigma=15$')

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.draw()

if GRAPH_DEBUG:
    plt.show() # wait until image is closed

# 2nd example

mu = 200
sigma = 25
x = mu + sigma*np.random.randn(10000)

fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(8, 4))

ax0.hist(x, 20, normed=1, histtype='stepfilled', facecolor='g', alpha=0.75)
ax0.set_title('stepfilled')

# Create a histogram by providing the bin edges (unequally spaced).
bins = [100, 150, 180, 195, 205, 220, 250, 300]
ax1.hist(x, bins, normed=1, histtype='bar', rwidth=0.8)
ax1.set_title('unequal bins')

plt.tight_layout()
plt.draw()

if GRAPH_DEBUG:
    plt.show() # wait until image is closed

# TODO: try to make plots with my numpy 1D histogram objects: is it possible, or do I have to use the "hist" method of matplotlib/pyplot?

# plot 2D histogram in 3D (with bars)
# ref: http://matplotlib.org/examples/mplot3d/hist3d_demo.html

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x, y = np.random.rand(2, 100) * 4
hist, xedges, yedges = np.histogram2d(x, y, bins=4)

elements = (len(xedges) - 1) * (len(yedges) - 1)
xpos, ypos = np.meshgrid(xedges[:-1]+0.25, yedges[:-1]+0.25)

xpos = xpos.flatten()
ypos = ypos.flatten()
zpos = np.zeros(elements)
dx = 0.5 * np.ones_like(zpos)
dy = dx.copy()
dz = hist.flatten()

ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_ylabel('HIST quantity')

plt.draw()

if GRAPH_DEBUG:
    plt.show() # wait until image is closed

# plot 2D histogram in 2D (with squares)
# not working!!! I guess I need to do a scatter plot!!!

#fig = plt.figure()
###ax = fig.add_subplot(111, projection='2d')
#ax = fig.add_subplot(111)
#x, y = np.random.rand(2, 100) * 4
#hist, xedges, yedges = np.histogram2d(x, y, bins=4)
#
##elements = (len(xedges) - 1) * (len(yedges) - 1)
#xpos, ypos = np.meshgrid(xedges[:-1]+0.25, yedges[:-1]+0.25)
#
#xpos = xpos.flatten()
#ypos = ypos.flatten()
#dx = hist.flatten()
#dy = dx.copy()
#
##ax.bar2d(xpos, ypos, dx, dy, color='b', zsort='average')
##ax.bar(xpos, ypos, dx, dy, color='b', zsort='average')
##ax.imshow([dx,dy], origin='lower')
#ax.set_xlabel('X')
#ax.set_ylabel('Y')
#
#plt.draw()
#
#if GRAPH_DEBUG:
#    plt.show() # wait until image is closed

# plot 2D histogram in 2D (scatter)
#(see a better example further below)

fig = plt.figure()
ax = fig.add_subplot(111)
x, y = np.random.rand(2, 100) * 4
hist, xedges, yedges = np.histogram2d(x, y, bins=4)

#elements = (len(xedges) - 1) * (len(yedges) - 1)
xpos, ypos = np.meshgrid(xedges[:-1]+0.25, yedges[:-1]+0.25)

xpos = xpos.flatten()
ypos = ypos.flatten()
dx = hist.flatten()
dy = dx.copy()
area = dx*dy

c = plt.scatter(xpos, ypos, c='b', s=area, cmap=plt.cm.hsv)
ax.set_xlabel('X')
ax.set_ylabel('Y')
plt.draw()

if GRAPH_DEBUG:
    plt.show() # wait until image is closed

# plot 2D histogram in 2D (with colors, similar colz in root)
# ref: http://matplotlib.org/examples/pylab_examples/hist2d_log_demo.html
# not exportable to 3D, because there is no hist3d!!!

fig = plt.figure()
ax = fig.add_subplot(111)

#normal distribution center at x=0 and y=5
x = randn(100000)
y = randn(100000)+5

plt.hist2d(x, y, bins=40, norm=LogNorm())
ax.set_xlabel('X')
ax.set_ylabel('Y')
plt.colorbar(label='represented quantity')
plt.draw()

if GRAPH_DEBUG:
    plt.show() # wait until image is closed

# For plotting 3D hist in box (cube) plot use:
#  x, y, z as middle points of bins defiend as edges (1 point less as bin edges) (bin centers)
#  data: cubes with side equal to the value of the histogram in the coresponding bin
#        eventually scaled by a constant factor for better plot visibility

# plot 3D histogram in 3D (with boxes)

data_z = np.array([21.,22.,23.,24.,25.]) # works for histogram2d
bins_3d = np.array([[0., 2., 4., 6.], [10., 12., 14., 16.], [20., 22., 24., 26.]])
astro_table['Z'] = data_z
print("astro_table")
print(astro_table)
print("structured data", data)
data = np.vstack([astro_table['X'], astro_table['Y'], astro_table['Z']]).T
print("non-structured data")
print(data)
hdd_table, edges_table = np.histogramdd(data, bins_3d)
print("edges_table")
print(repr(edges_table))
print("hdd_table")
print(repr(hdd_table))
   
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# define point positions as bin centers
x_edges_low = edges_table[0][:-1]
x_edges_high = edges_table[0][1:]
xpos = (x_edges_low + x_edges_high)/2.

y_edges_low = edges_table[1][:-1]
y_edges_high = edges_table[1][1:]
ypos = (y_edges_low + y_edges_high)/2.

z_edges_low = edges_table[2][:-1]
z_edges_high = edges_table[2][1:]
zpos = (z_edges_low + z_edges_high)/2.

# define grid of points (i.e. coords for each 2D bin)
xpos, ypos, zpos = np.meshgrid(xpos, ypos, zpos)
xpos = xpos.flatten() # reduce to a list of coords
ypos = ypos.flatten() # reduce to a list of coords
zpos = zpos.flatten() # reduce to a list of coords

# dimensions of the boxes
dx = hdd_table.flatten()
#dx *= 0.5 # scale! # WARNING! if I scale, the boxes are shifted w.r.t. the bin centers!!! TODO: FIXME!!!!
dy = dx.copy()
dz = dx.copy()

ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.draw()

if GRAPH_DEBUG:
    plt.show() # wait until image is closed

# For plotting 2D hist in scatter plot use:
#  x, y as middle points of bins defiend as edges (1 point less as bin edges) (bin centers)
#  data: circles with radius equal to the value of the histogram in the coresponding bin
#        eventually scaled by a constant factor for better plot visibility
# Analogous procedure for 3D hist scatter plot.

# plot 2D histogram in 2D (scatter)
# ref: http://matplotlib.org/examples/pie_and_polar_charts/polar_scatter_demo.html

fig = plt.figure()
ax = fig.add_subplot(111)

# define point positions as bin centers
x_edges_low = x_edges_data[:-1]
x_edges_high = x_edges_data[1:]
xpos = (x_edges_low + x_edges_high)/2.

y_edges_low = y_edges_data[:-1]
y_edges_high = y_edges_data[1:]
ypos = (y_edges_low + y_edges_high)/2.

# define grid of points (i.e. coords for each 2D bin)
xpos, ypos = np.meshgrid(xpos, ypos)
xpos = xpos.flatten() # reduce to a list of coords
ypos = ypos.flatten() # reduce to a list of coords

# area for the circles in the scatter plot
r = h2d_data.flatten() # reduce to a list of bin contents
area = r**2
area *= 100. # scale!

colormap = plt.cm.hsv
# ref: http://matplotlib.org/examples/color/colormaps_reference.html
# beware: scatter is called in "ax" for 3D plots, but in "plt" for 2D plots!
c = plt.scatter(xpos, ypos, c='b', s=area, cmap=colormap)
ax.set_xlabel('X')
ax.set_ylabel('Y')
plt.draw()

if GRAPH_DEBUG:
    plt.show() # wait until image is closed

# plot 3D histogram in 3D (scatter)
# ref: http://matplotlib.org/examples/mplot3d/scatter3d_demo.html

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# define point positions as bin centers
x_edges_low = edges_table[0][:-1]
x_edges_high = edges_table[0][1:]
xpos = (x_edges_low + x_edges_high)/2.

y_edges_low = edges_table[1][:-1]
y_edges_high = edges_table[1][1:]
ypos = (y_edges_low + y_edges_high)/2.

z_edges_low = edges_table[2][:-1]
z_edges_high = edges_table[2][1:]
zpos = (z_edges_low + z_edges_high)/2.

# define grid of points (i.e. coords for each 2D bin)
xpos, ypos, zpos = np.meshgrid(xpos, ypos, zpos)
xpos = xpos.flatten() # reduce to a list of coords
ypos = ypos.flatten() # reduce to a list of coords
zpos = zpos.flatten() # reduce to a list of coords

# area for the circles in the scatter plot
r = hdd_table.flatten() # reduce to a list of bin contents
area = r**2
area *= 100. # scale!

colormap = plt.cm.hsv
# ref: http://matplotlib.org/examples/color/colormaps_reference.html
# beware: scatter is called in "ax" for 3D plots, but in "plt" for 2D plots!
c = ax.scatter(xpos, ypos, zpos, c='b', s=area, cmap=colormap)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.draw()

if GRAPH_DEBUG:
    plt.show() # wait until image is closed

plt.show() #don't quit at the end

from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
from gammapy.background import CubeBackgroundModel

# test read table and image cube bg model files
group = 1
inputdir = 'bg_cube_models_new_write/'
infile = inputdir + 'bg_cube_model_group{}'.format(group)

# read table file
bg_cube_model_table = CubeBackgroundModel.read('{}_table.fits.gz'.format(infile), format='table')

# read image file
bg_cube_model_image = CubeBackgroundModel.read('{}_image.fits.gz'.format(infile), format='image')

assert (bg_cube_model_image.background_cube.data == bg_cube_model_table.background_cube.data).all()

import IPython; IPython.embed()

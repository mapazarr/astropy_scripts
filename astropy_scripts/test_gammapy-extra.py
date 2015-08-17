from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
from tempfile import NamedTemporaryFile
from numpy.testing import assert_allclose
from gammapy.background import Cube
from gammapy import datasets

CACHE = 1 # set to 0 to ignore downloaded files in the cache
          # (or use rm ~/.astropy/cache)

filename = '../test_datasets/background/bg_cube_model_test.fits'
filename = datasets.get_path(filename, location='remote', cache=CACHE)
bg_cube_model_1 = Cube.read(filename, format='table', scheme='bg_cube')

outfile = NamedTemporaryFile(suffix='.fits').name
bg_cube_model_1.write(outfile, format='image')

# test if values are correct in the saved file: compare both files
bg_cube_model_2 = Cube.read(outfile, format='image', scheme='bg_cube')
assert_allclose(bg_cube_model_2.data,
                bg_cube_model_1.data)
assert_allclose(bg_cube_model_2.coordx_edges,
                bg_cube_model_1.coordx_edges)
assert_allclose(bg_cube_model_2.coordy_edges,
                bg_cube_model_1.coordy_edges)
assert_allclose(bg_cube_model_2.energy_edges,
                bg_cube_model_1.energy_edges)

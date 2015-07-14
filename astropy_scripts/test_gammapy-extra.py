from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
from tempfile import NamedTemporaryFile
from numpy.testing import assert_allclose
from gammapy.background import CubeBackgroundModel
from gammapy import datasets

filename = '../test_datasets/background/bg_cube_model_test.fits'
filename = datasets.get_path(filename, location='remote')
bg_model_1 = CubeBackgroundModel.read(filename, format='table')

outfile = NamedTemporaryFile(suffix='.fits').name
bg_model_1.write(outfile, format='image')

# test if values are correct in the saved file: compare both files
bg_model_2 = CubeBackgroundModel.read(outfile, format='image')
assert_allclose(bg_model_2.background,
                bg_model_1.background)
assert_allclose(bg_model_2.detx_bins, # fails if ouside from gammapy!
                bg_model_1.detx_bins)
assert_allclose(bg_model_2.dety_bins, # fails if ouside from gammapy!
                bg_model_1.dety_bins)
assert_allclose(bg_model_2.energy_bins,
                bg_model_1.energy_bins)

# the error is produced because the wrong (outdated) version of the
# file is fetched from the gammapy-extra repo!

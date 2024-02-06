import copy
import glob
import netCDF4
import numpy as np
import os
import subprocess
import glob
import unittest

import xios_examples.shared_testing as xshared

this_path = os.path.realpath(__file__)
this_dir = os.path.dirname(this_path)

class TestResampleAxis(xshared._TestCase):
    test_dir = this_dir
    transient_inputs = ['axis_input.nc']
    transient_outputs = ['axis_output.nc']


# A list of input `.cdl` files where XIOS is known to produce different
# output from the expected output data
# for future investigation / ToDo
known_failures = ['test_axis_input_edge_simple_square_ten']

# iterate through `.cdl` files in this test case folder
for f in glob.glob('{}/*.cdl'.format(this_dir)):
    # unique name for the test
    tname = 'test_{}'.format(os.path.splitext(os.path.basename(f))[0])
    # add the test as an attribute (function) to the test class
    if tname in known_failures:
        # set decorator @unittest.expectedFailure
        setattr(TestResampleAxis, tname,
                unittest.expectedFailure(TestResampleAxis.make_a_resample_test(f)))
    else:
        setattr(TestResampleAxis, tname,
                TestResampleAxis.make_a_resample_test(f))

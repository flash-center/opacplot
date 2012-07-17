import os

import nose
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, assert_in, \
                       assert_true, with_setup

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from opacplot import ionmix, undo_ionmix

def parser_teardown():
    files = os.listdir('.')
    h5files = [f for f in files if f.endswith('.h5')]
    for f in h5files:
        os.remove(f)
    cn4files = [f for f in files if f.endswith('.cn4')]
    for f in cn4files:
        os.remove(f)

"""# @with_setup(lambda: None, parser_teardown)
def test_parse1():
    opp = ionmix.parse('imxtest.cn4')
"""

# @with_setup(lambda: None, parser_teardown)
def test_parse2():
    oppundo = undo_ionmix.parse('imxtest.h5')

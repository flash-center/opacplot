import os

import tables as tb
import numpy as np

BASIC_FILTERS = tb.Filters(complevel=5, complib='zlib', shuffle=True, fletcher32=False)

class OppFile(object):
    """The opacplot file."""

    def __init__(self, filename, format=None, *args):
        if os.path.isfile(filename) and not tb.isHDF5File(filename):
            # FIXME to be based on mode
            pass
        self._handle = tb.openFile(filename, mode='a', filters=BASIC_FILTERS)
        self.root = self._handle.root

    def __del__(self):
        self._handle.close()

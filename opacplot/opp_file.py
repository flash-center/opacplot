import os

import tables as tb
import numpy as np

from .ionmix import parse
from .utils import isopacplot

BASIC_FILTERS = tb.Filters(complevel=5, complib='zlib', shuffle=True, fletcher32=False)

class OppFile(object):
    """The opacplot file."""

    def __init__(self, filename, format=None, *args, **kwargs):
        if not isopacplot(filename, format):
            opp = parse(filename, *args, **kwargs)
            self.__dict__.update(opp.__dict__)
            del opp
            return 
        self._handle = tb.openFile(filename, mode='a', filters=BASIC_FILTERS)
        self.root = self._handle.root
        self.__file__ = filename

    def __del__(self):
        self._handle.close()

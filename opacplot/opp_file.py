import os
import sys

import tables as tb
import numpy as np

from .utils import isopacplot

BASIC_FILTERS = tb.Filters(complevel=5, complib='zlib', shuffle=True, fletcher32=False)

FORMAT_PARSERS = {
    'ionmix': ('ionmix', 'parse'),
    }

class OppFile(object):
    """The opacplot file."""

    def __init__(self, filename, format=None, *args, **kwargs):
        if not isopacplot(filename, format):
            opp = self._parse(filename, format, *args, **kwargs)
            self.__dict__.update(opp.__dict__)
            del opp
            return 
        self._handle = tb.openFile(filename, mode='a', filters=BASIC_FILTERS)
        self.root = self._handle.root
        self.__file__ = filename

    @classmethod
    def _parse(cls, filename, format, *args, **kwargs):
        if format in FORMAT_PARSERS:
            opp = cls._format_parse(filename, format, *args, **kwargs)
        else:
            for form in FORMAT_PARSERS.keys():
                try:
                    opp = cls._format_parse(filename, form, *args, **kwargs)
                    break
                except:
                    continue
            else:
                msg = "format of {0} could not be determined from available types ({1})"
                msg = msg.format(filename, ", ".join(FORMAT_PARSERS.keys()))
                raise RuntimeError(msg)
        return opp

    @classmethod
    def _format_parse(cls, filename, format, *args, **kwargs):
        if '.' not in sys.path:
            sys.path.insert(0, '.')
        formatmod, parserfunc = FORMAT_PARSERS[format]
        formatmod = __import__(formatmod, globals(), locals(), fromlist=[None])
        parserfunc = getattr(formatmod, parserfunc)
        opp = parserfunc(filename, *args, **kwargs)
        return opp

    def __del__(self):
        self._handle.close()

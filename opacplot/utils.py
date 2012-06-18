import os

import tables as tb

def munge_h5filename(filename, h5filename):
    """Gets an opacplot hdf5 filename from existing names."""
    if not isinstance(h5filename, basestring):
        if isinstance(filename, basestring): 
            h5filename = filename.rpartition(".")[0] + '.h5'
        else:
            h5filename = "opp.h5"
    return h5filename


def isopacplot(filename, format=None):
    """Determines if a file is in of opacplot format."""
    return (format in [None, "opp"]) and isinstance(filename, basestring) and \
            ((os.path.isfile(filename) and tb.isHDF5File(filename)) or \
              not os.path.exists(filename))



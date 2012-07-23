from StringIO import StringIO
import numpy
import re
import math
import tables as tb
import itertools as it

from opl_grid import OplGrid
from constants import ERG_TO_JOULE

from ionmix import writeIonmixFile
from opp_file import OppFile
from utils import munge_cn4filename

# pick filters from pytables                                                                         
BASIC_FILTERS = tb.Filters(complevel=9, complib='zlib', shuffle=True, fletcher32=False)

def parse(h5filename, cn4filename=None, h5groupname=None, mpi=1.0, twot=False, man=True, verbose=False, *args, **kwargs):
    """Does some things like calling function to read hdf5 files and write cn4 files"""
    ionmixdata = _parse(h5filename, mpi, twot, man, verbose, *args, **kwargs)
    oppundo = writeIonmixFile(ionmixdata['filename'], ionmixdata['zvals'], ionmixdata['fracs'],
                              ionmixdata['ndens'], ionmixdata['ntemps'], ionmixdata['numdens'], 
                              ionmixdata['temps'], ionmixdata['zbar'], ionmixdata['dzdt'], 
                              ionmixdata['pion'], ionmixdata['pele'], ionmixdata['dpidt'],
                              ionmixdata['dpedt'], ionmixdata['eion'], ionmixdata['eele'],
                              ionmixdata['cvion'], ionmixdata['cvele'], ionmixdata['deidn'],
                              ionmixdata['deedn'], ionmixdata['ngroups'], ionmixdata['opac_bounds'], 
                              ionmixdata['rosseland'], ionmixdata['planck_absorb'], ionmixdata['planck_emiss'])

    return oppundo

def _parse(h5filename, mpi=1.0, cn4filename=None, twot=False, man=True, verbose=False, *args, **kwargs):
    """Read opp file, commits to memory"""
    if verbose:
        print "Reading h5 file '{0}'/n".format(h5filename)

    h5file = tb.openFile(h5filename, mode='a', filters=BASIC_FILTERS)

    ionmixdata = {}
    ionmixdata = {'filename': munge_cn4filename(h5filename, cn4filename)}

    names = [arr.name for arr in h5file.walkNodes('/Ionmix', classname='CArray')]
    namegroups = {'norm':[],'_dens':[],'_temps':[]}

    def endsorter(name):
        if name.endswith('_dens'):
            return "_dens"
        elif name.endswith('_temps'):
            return "_temps"
        else:
            return "norm"

    for key, group in it.groupby(names, endsorter):
        namegroups[key].extend(group)

    for info_type in h5file.walkNodes('/Ionmix', classname='CArray'):
        ionmixdata[info_type.name] = info_type.read()

    for info_type in h5file.walkNodes('/Ionmix', classname='EArray'):
        ionmixdata[info_type.name] = info_type.read()

    ionmixdata['temps'] = h5file.root.Ionmix.pion_temps
    ionmixdata['dens'] = h5file.root.Ionmix.pion_dens
    ionmixdata['ntemps'] = len(h5file.root.Ionmix.pion_temps)
    ionmixdata['ndens'] = len(h5file.root.Ionmix.pion_dens)
    ionmixdata['ngroups'] = h5file.root.Ionmix._f_getAttr('ngroups')
    ionmixdata['fracs'] = h5file.root.Ionmix.fracs    
    ionmixdata['zvals'] = h5file.root.Ionmix.zvals
    ionmixdata['opac_bounds'] = h5file.root.Ionmix.opac_bounds

    # check that "_temps" is list of files containing identical info
    for tempck in namegroups['_temps']:
        if all(ionmixdata[tempck]) != all(ionmixdata[namegroups['_temps'][0]]):
            msg = "temperature files {0} {1} not identical".format(ionmixdata[namegroups['_temps'][tempck]], ionmixdata[namegroups['_temps'][0]])
            raise RuntimeError(msg)

    # check that "_temps" is list of files containing identical info
    for densck in namegroups['_dens']:
        if all(ionmixdata[densck]) != all(ionmixdata[namegroups['_dens'][0]]):
            msg = "density files {0} {1} not identical".format(ionmixdata[namegroups['_dens'][densck]], ionmixdata[namegroups['_dens'][0]])
            raise RuntimeError(msg)

    return ionmixdata

from StringIO import StringIO
import numpy
import re
import math
import tables as tb
import itertools as it

from opl_grid import OplGrid
from constants import ERG_TO_JOULE

from opp_file import OppFile
from utils import munge_cn4filename

# pick filters from pytables                                                                         
BASIC_FILTERS = tb.Filters(complevel=9, complib='zlib', shuffle=True, fletcher32=False)

def parse(h5filename, cn4filename=None, h5groupname=None, mpi=1.0, twot=False, man=True, verbose=False, *args, **kwargs):
    """Does some things like calling function to read hdf5 files and write cn4 files"""
    ionmixdata = _parse(h5filename, mpi, twot, man, verbose, *args, **kwargs)
    oppundo = _tocn4(h5filename, ionmixdata, cn4filename, h5groupname)
    return oppundo

def _parse(h5filename, mpi=1.0, twot=False, man=True, verbose=False, *args, **kwargs):
    """Read opp file, commits to memory"""
    if verbose:
        print "Reading h5 file '{0}'/n".format(h5filename)

    h5file = tb.openFile(h5filename, mode='a', filters=BASIC_FILTERS)

    ionmixdata = {}

    for info_type in h5file.walkNodes('/imxtest', classname='CArray'):
        ionmixdata[info_type] = info_type

    names = [arr.name for arr in h5file.walkNodes('/imxtest', classname='CArray')]
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

    # check that "_temps" is list of files containing identical info
#    for tempck in namegroups['_temps']:
#        if ionmixdata[namegroups['_temps'][tempck]] != ionmixdata[namegroups['_temps'][1]]:
#            msg = "temperature files {0} {1} not identical".format(ionmixdata[namegroups['_temps'][tempck]], ionmixdata[namegroups['_temps'][1]])
#            raise RuntimeError(msg)

    # check that "_temps" is list of files containing identical info
#    for densck in namegroups['_dens']:
#        if ionmixdata[namegroups['_dens'][densck]] != ionmixdata[namegroups['_dens'][1]]:
#            msg = "density files {0} {1} not identical".format(ionmixdata[namegroups['_dens'][densck]], ionmixdata[namegroups['_dens'][1]])
#            raise RuntimeError(msg)

    ionmixdata['ntemp'] = len(h5file.root.imxtest.etot_temps)
    ionmixdata['ndens'] = len(h5file.root.imxtest.etot_dens)
    ionmixdata['ngroups'] = h5file.root.imxtest._f_getAttr('ngroups')
    ionmixdata['relative_fractions'] = h5file.root.imxtest.relative_fractions
    ionmixdata['atom_nums'] = h5file.root.imxtest.atomic_numbers
    ionmixdata['energy_group_bounds'] = h5file.root.imxtest.energy_group_bounds

    return ionmixdata



def _tocn4(h5filename, ionmixdata, cn4filename=None, twot=False):
#    if isinstance(h5filename, basestring):
#        oppun = OppFile(h5filename)
#    elif 

    if 'dzdt' in ionmixdata:
        twot = True

    cn4filename = munge_cn4filename(h5filename,cn4filename)

    cn4f = open(cn4filename, 'w')

    cn4f.write('          ')
    cn4f.write(str(ionmixdata['ntemp']))
    cn4f.write('          ')
    cn4f.write(str(ionmixdata['ndens']))
    cn4f.write('\n')
    cn4f.write(' atomic #s of gases:')
    cn4f.write('         ')
    cn4f.write(str(ionmixdata['atom_nums'][0]).strip('[]'))
    cn4f.write('         ')
    cn4f.write(str(ionmixdata['atom_nums'][1]).strip('[]'))
    cn4f.write('\n')
    cn4f.write(' relative fractions:')
    cn4f.write('  ')
    cn4f.write('{0:0.02e}'.format(float(str(ionmixdata['relative_fractions'][0]).strip('[]')),))
    cn4f.write('  ')
    cn4f.write('{0:0.02e}'.format(float(str(ionmixdata['relative_fractions'][1]).strip('[]')),))
    cn4f.write('\n')
    cn4f.write('         ')
    cn4f.write(str(ionmixdata['ngroups']))
    cn4f.write('\n')

    if twot:
        cn4f.write(ionmixdata['dzdt'])
        cn4f.write(ionmixdata['pion'])
        cn4f.write(ionmixdata['pele'])
        cn4f.write(ionmixdata['dpidt'])
        cn4f.write(ionmixdata['dpedt'])
        cn4f.write(ionmixdata['eion'])
        cn4f.write(ionmixdata['eele'])
        cn4f.write(ionmixdata['cvion'])
        cn4f.write(ionmixdata['cvele'])
        cn4f.write(ionmixdata['deidn'])
        cn4f.write(ionmixdata['deedn'])
        cn4f.write(ionmixdata['energy_group_bounds'])
        cn4f.write(ionmixdata['rosseland'])
        cn4f.write(ionmixdata['planck_absorb'])
        cn4f.write(ionmixdata['planck_emiss'])
    else:
        cn4f.write(str(ionmixdata['etot']))
        cn4f.write(ionmixdata['cvtot'])
        cn4f.write(ionmixdata['dedn'])
        cn4f.write(ionmixdata['energy_group_bounds'])
        cn4f.write(ionmixdata['rosseland'])
        cn4f.write(ionmixdata['planck_absorb'])
        cn4f.write(ionmixdata['planck_emiss'])

    cn4f.close()

#_tocn4('hello.h5')
#_tocn4(opp)

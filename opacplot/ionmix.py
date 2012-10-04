from StringIO import StringIO
import numpy as np
import re 
import math
import tables as tb

from opl_grid import OplGrid
from constants import ERG_TO_JOULE

from opp_file import OppFile
from utils import munge_h5filename, munge_h5groupname

# pick filters from pytables
# BASIC_FILTERS = tb.Filters(complevel=5, complib='zlib', shuffle=True, fletcher32=False)

joules_to_ergs = 1.0e+07

def parse(filename, h5filename='eos_opac_lib.h5', h5groupname=None, mpi=1.0, twot=True, man=True, verbose=False, *args, **kwargs):
    """Calls _parse to read in data from ionmix data file and commit to memory in dictionary. Thens calls _tohdf5 to generate or open existing h5 file, make modifications, exit."""
    ionmixdata = _parse(filename, mpi, twot, man, verbose, *args, **kwargs)
    opp = _tohdf5(filename, h5filename, h5groupname, ionmixdata)
    return opp

def _parse(filename, mpi=1.0, twot=True, man=True, verbose=False, *args, **kwargs):    
    ionmixdata = {}
    
    if verbose:
        print "Reading IONMIX file '{0}'\n".format(filename)

    f = open(filename, 'r')

    # Read the number of temperatures/densities:
    ntemp = int(f.read(10))
    ndens = int(f.read(10))
    ionmixdata['ntemp'] = ntemp
    ionmixdata['ndens'] = ndens
    
    # skip the rest of the line
    f.readline()

    # Read atomic numbers and mass fractions:
    ionmixdata['atomic_nums']=()
    f.read(21)
    fr = f.readline()
    for an in fr.split():
        ionmixdata['atomic_nums'] += (an,)

    ionmixdata['material_fracs']=()
    f.read(21)
    fr = f.readline()
    for rf in fr.split():
        ionmixdata['material_fracs'] += (rf,)

    # Setup temperature/density grid:
    if not man:
        # Read information about the temperature/density grid:
        ddens_log10 = float(f.read(12))
        dens0_log10 = float(f.read(12))
        dtemp_log10 = float(f.read(12))
        temp0_log10 = float(f.read(12))

        # Compute number densities:
        num_dens = np.logspace(dens0_log10, 
                               dens0_log10 + ddens_log10 * (ndens - 1), 
                               ndens)

        temps = np.logspace(temp0_log10, 
                            temp0_log10 + dtemp_log10 * (ntemp - 1), 
                            ntemp)

        # Read number of groups:
        ngroups = int(f.read(12))
    else:
        ngroups = int(f.read(12))

    # Read the rest of the file, remove all of the white space,
    # and store the string in self.data:
    data = StringIO(re.sub(r'\s', '', f.read()))
    f.close()

    if man:
        # For files where temperatures/densities are manually
        # specified, read the manual values here.
        temps = get_block(data, ntemp)
        num_dens = get_block(data, ndens)

    dens = num_dens * mpi

    if verbose: 
        print "  Number of temperatures: {0}".format(ntemp)
        for i in range(ntemp):
            print "%6i%27.16e" % (i, temps[i])

        print "\n  Number of densities: {0}".format(ndens)
        for i in range(ndens):
            print "%6i%21.12e%27.16e" % (i, dens[i], num_dens[i])
 
    ionmixdata['num_dens'] = num_dens
    ionmixdata['ngroups'] = ngroups
    ionmixdata['temps'] = temps
    ionmixdata['dens'] = dens
    ionmixdata['eos'] = read_eos(data, twot, ntemp, ndens, ngroups)
    ionmixdata['opac'] = read_opac(data, ntemp, ndens, ngroups, verbose)

    return ionmixdata


def get_block(data, n):
    arr = np.zeros(n)
    for i in range(n):
        arr[i] = float(data.read(12))
    return arr

def read_eos(data, twot, nt, nd, ng):
    eos = {}
    eos['avg_ionization'] = get_block(data, nd*nt).reshape(nd,nt)     # average charge state

    if twot:
        # Read in pressure, specific internal energies and
        # specific heats, but convert from J to ergs:
        eos['dcharge_dtemp']  = get_block(data, nd*nt).reshape(nd,nt)
        eos['ion_press']  = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
        eos['ele_press']  = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
        eos['dion_press_dion_temp'] = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
        eos['dele_press_dele_temp'] = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
        eos['ion_energy'] = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
        eos['ele_energy'] = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
        eos['dion_energy_dion_temp'] = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
        eos['dele_energy_dele_temp'] = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
        eos['dion_energy_dion_dens'] = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
        eos['dele_energy_dele_dens'] = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
            
    else: 
        # Read in e and cv, but convert from J to ergs:
        eos['etot']  = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs   #
        eos['cvtot'] = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs   #
        eos['dedn']  = get_block(data, nd*nt).reshape(nd,nt)   #

    return eos

def read_opac(data, nt, nd, ng, verbose):
    """
    Load the opacities from the file. The opacities are arranged
    in the file so that temperature varies the fastest, then
    density, and group number varies the slowest. Note, that this
    is not the ordering of the arrays once they are loaded.
    """
    
    opac = {}

    # Read group bounds in eV and convert to ergs:
    opac_bounds = get_block(data, ng+1)

    if verbose:
        print "\n  Number of Energy Groups: {0}".format(ng)
        for i in range(0, ng+1):
            print "%6i%15.6e" % (i, opac_bounds[i])

    rosseland     = np.empty((nd,nt,ng))
    planck_absorb = np.empty((nd,nt,ng))
    planck_emiss  = np.empty((nd,nt,ng))
        
    arr_ro = get_block(data, nd*nt*ng)
    arr_pa = get_block(data, nd*nt*ng)
    arr_pe = get_block(data, nd*nt*ng)

    i = 0
    for g in xrange(ng):
        for d in xrange(nd):
            for t in xrange(nt):
                rosseland[d,t,g]     = arr_ro[i]
                planck_absorb[d,t,g] = arr_pa[i]
                planck_emiss[d,t,g]  = arr_pe[i]
                i += 1    
                     
    opac['energy_group_bnds'] = opac_bounds
    opac['rosseland'] = rosseland
    opac['planck_absorb'] = planck_absorb
    opac['planck_emiss'] = planck_emiss
     
    return opac


def _tohdf5(filename, h5filename, h5groupname, ionmixdata):
    if h5filename != 'eos_opac_lib.h5':
        h5filename = munge_h5filename(filename, h5filename)
    
    opp = OppFile(h5filename)
    opph = opp._handle

    h5groupname = filename
    ionmix_group = opph.createGroup(opp.root, h5groupname)

    opph.setNodeAttr(ionmix_group, 'ngroups', ionmixdata['ngroups'], name=None)
    opph.setNodeAttr(ionmix_group, 'code', 'Ionmix', name=None) 
    opph.setNodeAttr(ionmix_group, 'ntemp', ionmixdata['ntemp'], name=None)
    opph.setNodeAttr(ionmix_group, 'ndens', ionmixdata['ndens'], name=None)




    ndset = opph.createCArray(ionmix_group, "num_dens", tb.Float64Atom(), ionmixdata['num_dens'].shape)
    for x in range(len(ionmixdata['num_dens'])):
        ndset[x] = ionmixdata['num_dens'][x]

    zvset = opph.createCArray(ionmix_group, "atomic_nums", tb.Int16Atom(), (len(ionmixdata['atomic_nums']), 1))
    for x in range(len(ionmixdata['atomic_nums'])):
        zvset[x] = ionmixdata['atomic_nums'][x]

    rfset = opph.createCArray(ionmix_group, "material_fracs", tb.Float64Atom(), (len(ionmixdata['material_fracs']), 1))
    for x in range(len(ionmixdata['material_fracs'])):
        rfset[x] = ionmixdata['material_fracs'][x]

    for key, value in ionmixdata['eos'].items():
        dset = opph.createCArray(ionmix_group, key, tb.Float64Atom(), ionmixdata['eos'][key].shape)
        dset[:] = value
        dset.attrs.dets = ("dens", "temps")
        set_dens = opph.createCArray(ionmix_group, key + "_dens", tb.Float64Atom(), ionmixdata['dens'].shape)
        set_dens[:] = ionmixdata['dens']
        set_temps = opph.createCArray(ionmix_group, key + "_temps", tb.Float64Atom(), ionmixdata['temps'].shape)
        set_temps[:] = ionmixdata['temps']
        
    for key, value in ionmixdata['opac'].items():
        dset = opph.createCArray(ionmix_group, key, tb.Float64Atom(), ionmixdata['opac'][key].shape)
        dset[:] = value
        if key != 'opac_bounds':
            dset.attrs.dets = ("dens", "temps")
            set_dens = opph.createCArray(ionmix_group, key + "_dens", tb.Float64Atom(), ionmixdata['dens'].shape)
            set_dens[:] = ionmixdata['dens']
            set_temps = opph.createCArray(ionmix_group, key + "_temps", tb.Float64Atom(), ionmixdata['temps'].shape)
            set_temps[:] = ionmixdata['temps']

    h = open(filename, 'r')
    content = h.read()
    info = opph.createGroup(ionmix_group, "info")

    a = tb.StringAtom(itemsize=1)
    # Use 'a' as the object type for the enlargeab3le array
    origfile = opph.createEArray(info, filename, a, (0,), "Chars")
    origfile.append([c for c in content])

    opph.flush()

    return opp


def load_h5file(h5filename, cn4filename=None, h5groupname=None, mpi=1.0, twot=False, man=True, verbose=False, *args, **\
kwargs):
    """Does some things like calling function to read hdf5 files and write cn4 files"""
    ionmixdata = _load_h5file(h5filename, mpi, twot, man, verbose, *args, **kwargs)
    oppundo = write_ionmix_file(ionmixdata['filename'], ionmixdata['atomic_nums'], ionmixdata['material_fracs'],
                                ionmixdata['ndens'], ionmixdata['ntemp'], ionmixdata['num_dens'],
                                ionmixdata['temps'], ionmixdata['avg_ionization'],
                                ionmixdata['dcharge_dtemp'],
                                ionmixdata['ion_press'], ionmixdata['ele_press'],
                                ionmixdata['dion_press_dion_temp'], ionmixdata['dele_press_dele_temp'],
                                ionmixdata['ion_energy'], ionmixdata['ele_energy'],
                                ionmixdata['dion_energy_dion_temp'], ionmixdata['dele_energy_dele_temp'],
                                ionmixdata['dion_energy_dion_dens'], ionmixdata['dele_energy_dele_dens'],
                                ionmixdata['ngroups'], ionmixdata['energy_group_bnds'],
                                ionmixdata['rosseland'], ionmixdata['planck_absorb'],
                                ionmixdata['planck_emiss'])

    return oppundo

def _load_h5file(h5filename, mpi=1.0, cn4filename=None, twot=False, man=True, verbose=False, *args, **kwargs):
    """Read opp file, commits to memory"""
    if verbose:
        print "Reading h5 file '{0}'/n".format(h5filename)

    h5file = tb.openFile(h5filename, mode='r')

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
    ionmixdata['ntemp'] = len(h5file.root.Ionmix.pion_temps)
    ionmixdata['ndens'] = len(h5file.root.Ionmix.pion_dens)
    ionmixdata['ngroups'] = h5file.root.Ionmix._f_getAttr('ngroups')
    ionmixdata['material_fracs'] = h5file.root.Ionmix.material_fracs
    ionmixdata['atomic_nums'] = h5file.root.Ionmix.atomic_nums
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



def write_ionmix_file(fn, atomic_nums, material_fracs, ndens, ntemp, num_dens, temps, 
                      zbar=None,  dzdt=None, pion=None, pele=None,
                      dpidt=None, dpedt=None, eion=None, eele=None,
                      cvion=None, cvele=None, deidn=None, deedn=None,
                      ngroups=None, opac_bounds=None,
                      rosseland=None, planck_absorb=None, planck_emiss=None):

    if avg_ionization == None:  zbar = np.zeros((ndens,ntemp))
    if dcharge_dtemp == None:  dzdt = np.zeros((ndens,ntemp))
    if ion_press == None:  pion = np.zeros((ndens,ntemp))
    if ele_press == None:  pele = np.zeros((ndens,ntemp))
    if dion_press_dion_temp == None: dpidt = np.zeros((ndens,ntemp))
    if dele_press_dele_temp == None: dpedt = np.zeros((ndens,ntemp))
    if ion_energy == None:  eion = np.zeros((ndens,ntemp))
    if ele_energy == None:  eele = np.zeros((ndens,ntemp))
    if dion_energy_dion_temp == None: cvion = np.zeros((ndens,ntemp))
    if dele_energy_dion_dens == None: cvele = np.zeros((ndens,ntemp))
    if dion_energy_dion_dens == None: deidn = np.zeros((ndens,ntemp))
    if dele_energy_dele_dens == None: deedn = np.zeros((ndens,ntemp))

    if ngroups       == None: ngroups = 1
    if energy_group_bnds   == None: opac_bounds = (0.0,1.0)
    if rosseland     == None: rosseland = np.zeros((ndens,ntemp,ngroups))
    if planck_absorb == None: planck_absorb = np.zeros((ndens,ntemp,ngroups))
    if planck_emiss  == None: planck_emiss = np.zeros((ndens,ntemp,ngroups))

    # Write the header:
    f = open(fn,'w')
    f.write("%10i%10i\n" % (ntemp,ndens))
    f.write(" atomic #s of gases: ")
    for z in atomic_nums: f.write("%10i" % z)
    f.write("\n relative fractions: ")
    for frac in material_fracs: f.write("%10.2E" % frac)
    f.write("\n")
    f.write("%12i\n" % (ngroups))

    # Write temperature/density grid and number of groups:
    def convert(num):
        string_org = "%12.5E" % (num)
        negative = (string_org[0] == "-")            
        lead = "-." if negative else "0."
        string = lead + string_org[1] + string_org[3:8] + "E"
    
        # Deal with the exponent:
        
        # Check for zero:
        if int(string_org[1] + string_org[3:8]) == 0:
            return string + "+00"

        # Not zero:
        expo = int(string_org[9:]) + 1
        if expo < 0:
            string += "-"
        else:
            string += "+"
        string += "%02d" % abs(expo)
        return string

    def write_block(var):
        count = 0
        for n in xrange(len(var)):
            count += 1

            f.write("%s" % convert(var[n]))
                
            if count == 4:
                count = 0
                f.write("\n")

        if count != 0: f.write("\n")

    def write_opac_block(var):
        count = 0
        for g in xrange(ngroups):
            for jd in xrange(ndens):
                for jt in xrange(ntemp):
                    count += 1

                    f.write("%s" % convert(var[jd,jt,g]))

                    if count == 4:
                        count = 0
                        f.write("\n")

        if count != 0: f.write("\n")

    write_block(temps)
    write_block(num_dens)

    write_block(avg_ionization.flatten())    
    write_block(dcharge_dtemp.flatten())
    write_block(ion_press.flatten()*ERG_TO_JOULE)
    write_block(ele_press.flatten()*ERG_TO_JOULE)
    write_block(dion_press_dion_temp.flatten()*ERG_TO_JOULE)
    write_block(dele_press_dele_temp.flatten()*ERG_TO_JOULE)
    write_block(ion_energy.flatten()*ERG_TO_JOULE)
    write_block(ele_energy.flatten()*ERG_TO_JOULE)
    write_block(dion_energy_dion_temp.flatten()*ERG_TO_JOULE)
    write_block(dele_energy_dele_temp.flatten()*ERG_TO_JOULE)
    write_block(dion_energy_dion_dens.flatten()*ERG_TO_JOULE)
    write_block(dele_energy_dele_dens.flatten()*ERG_TO_JOULE)

    write_block(energy_group_bnds)
    write_opac_block(rosseland)
    write_opac_block(planck_absorb)
    write_opac_block(planck_emiss)


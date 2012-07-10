from StringIO import StringIO
import numpy as np
import re 
import math

from opl_grid import OplGrid
from constants import ERG_TO_JOULE

from opp_file import OppFile
from utils import munge_h5filename

joules_to_ergs = 1.0e+07

def parse(filename, h5filename=None, mpi=0.0, twot=False, man=True, verbose=False, *args, **kwargs):
    if verbose:
        print "Reading IONMIX file '{0}'\n".format(filename)

    f = open(filename, 'r')

    # Read the number of temperatures/densities:
    ntemp = int(f.read(10))
    ndens = int(f.read(10))

    # Skip the next three lines:
    for i in range(3):
        f.readline()

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
        f.readline()

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

    eos = read_eos(data, twot, ntemp, ndens, ngroups)
    opac = read_opac(data, ntemp, ndens, ngroups, verbose)

    h5filename = munge_h5filename(filename, h5filename)
    opp = OppFile(h5filename)
    return opp


def get_block(data, n):
    arr = np.zeros(n)
    for i in range(n):
        arr[i] = float(data.read(12))
    return arr

def read_eos(data, twot, nt, nd, ng):

    eos = {}
    
    eos['zbar'] = get_block(data, nd*nt).reshape(nd,nt)

    if twot:
        # Read in pressure, specific internal energies and
        # specific heats, but convert from J to ergs:
        eos['dzdt']  = get_block(data, nd*nt).reshape(nd,nt)
        eos['pion']  = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
        eos['pele']  = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
        eos['dpidt'] = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
        eos['dpedt'] = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
        eos['eion']  = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
        eos['eele']  = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
        eos['cvion'] = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
        eos['cvele'] = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
        eos['deidn'] = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
        eos['deedn'] = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
            
    else: 
        # Read in e and cv, but convert from J to ergs:
        eos['etot']  = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
        eos['cvtot'] = get_block(data, nd*nt).reshape(nd,nt) * joules_to_ergs
        eos['dedn']  = get_block(data, nd*nt).reshape(nd,nt)

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
                     
    opac['rosseland'] = rosseland
    opac['planck_absorb'] = planck_absorb
    opac['planck_emiss'] = planck_emiss
     
    return opac
        
def oplAbsorb(dens, temps, opac_bounds):
    return OplGrid(dens, temps, opac_bounds, 
                   lambda jd, jt: planck_absorb[jd,jt,:])

def oplEmiss(dens, temps, opac_bounds):
    return OplGrid(dens, temps, opac_bounds, 
                   lambda jd, jt: planck_emiss[jd,jt,:])

def oplRosseland(dens, temps, opac):
    return OplGrid(dens, temps, opac_bounds, 
                   lambda jd, jt: rosseland[jd,jt,:])

def write(self, fn, zvals, fracs, twot=None, man=None):
    if twot == None: twot = twot
    if twot == True and twot == False:
        raise ValueError("Error: Cannot write two-temperature data")

    if man == None: man = man
    if man == False and man == True:
        raise ValueError("Error: Cannot write manual temp/dens points")

    # Write the header:
    f = open(fn,'w')
    f.write("%10i%10i\n" % (ntemp, ndens))
    f.write(" atomic #s of gases: ")
    for z in zvals: f.write("%10i" % z)
    f.write("\n relative fractions: ")
    for frac in fracs: f.write("%10.2E" % frac)
    f.write("\n")

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

    if not man:    
        f.write("%s%s%s%s" % (convert(ddens_log10), 
                              convert(dens0_log10), 
                              convert(dtemp_log10), 
                              convert(temp0_log10)) )
            
    f.write("%12i\n" % ngroups)

    if man == True:
        write_block(temps)
        write_block(numDens)

    write_block(zbar.flatten())

    if twot == False:
        write_block(etot.flatten()/joules_to_ergs)
        write_block(cvtot.flatten()/joules_to_ergs)
        write_block(enntab.flatten())

    else:
        write_block(dzdt.flatten())
        write_block(pion.flatten()/joules_to_ergs)
        write_block(pele.flatten()/joules_to_ergs)
        write_block(dpidt.flatten()/joules_to_ergs)
        write_block(dpedt.flatten()/joules_to_ergs)
        write_block(eion.flatten()/joules_to_ergs)
        write_block(eele.flatten()/joules_to_ergs)
        write_block(cvion.flatten()/joules_to_ergs)
        write_block(cvele.flatten()/joules_to_ergs)
        write_block(deidn.flatten()/joules_to_ergs)
        write_block(deedn.flatten()/joules_to_ergs)

    write_block(opac_bounds)
    write_opac_block(rosseland)
    write_opac_block(planck_absorb)
    write_opac_block(planck_emiss)


def extendToZero(self, nt, nd, ng):
    """
    This routine adds another temperature point at zero
    """
    
    arr = np.zeros((nd,nt+1))
    arr[:,1:] = dzdt[:,:]
    dzdt = arr

    arr = np.zeros((nd,nt+1))
    arr[:,1:] = pion[:,:]
    pion = arr

    arr = np.zeros((nd,nt+1))
    arr[:,1:] = pele[:,:]
    pele = arr

    arr = np.zeros((nd,nt+1))
    arr[:,1:] = dpidt[:,:]
    dpidt = arr

    arr = np.zeros((nd,nt+1))
    arr[:,1:] = dpedt[:,:]
    dpedi = arr

    arr = np.zeros((nd,nt+1))
    arr[:,1:] = eion[:,:]
    eion = arr

    arr = np.zeros((nd,nt+1))
    arr[:,1:] = eele[:,:]
    eele = arr

    arr = np.zeros((nd,nt+1))
    arr[:,1:] = zbar[:,:]
    zbar = arr

    arr = np.zeros((nd,nt+1))
    arr[:,1:] = cvion[:,:]
    cvion = arr

    arr = np.zeros((nd,nt+1))
    arr[:,1:] = self.cvele[:,:]
    self.cvele = arr

    arr = np.zeros((nd,nt+1))
    arr[:,1:] = deidn[:,:]
    deidn = arr

    arr = np.zeros((nd,nt+1))
    arr[:,1:] = deedn[:,:]
    deedn = arr

    arr = np.zeros((nd,nt+1,ng))
    arr[:,1:,:] = rosseland[:,:,:]
    rosseland = arr

    arr = np.zeros((nd,nt+1,ng))
    arr[:,1:,:] = planck_absorb[:,:,:]
    planck_absorb = arr

    arr = np.zeros((nd,nt+1,ng))
    arr[:,1:,:] = planck_emiss[:,:,:]
    planck_emiss = arr

    # Reset temperatures:

    arr = np.zeros((nt+1))
    arr[1:] = temps[:]
    temps = arr

    ntemp += 1


def writeIonmixFile(fn, zvals, fracs, ndens, ntemps, numDens, temps, 
                    zbar=None,  dzdt=None, pion=None, pele=None,
                    dpidt=None, dpedt=None, eion=None, eele=None,
                    cvion=None, cvele=None, deidn=None, deedn=None,
                    ngroups=None, opac_bounds=None,
                    rosseland=None, planck_absorb=None, planck_emiss=None):

    if  zbar == None:  zbar = np.zeros((ndens,ntemps))
    if  dzdt == None:  dzdt = np.zeros((ndens,ntemps))
    if  pion == None:  pion = np.zeros((ndens,ntemps))
    if  pele == None:  pele = np.zeros((ndens,ntemps))
    if dpidt == None: dpidt = np.zeros((ndens,ntemps))
    if dpedt == None: dpedt = np.zeros((ndens,ntemps))
    if  eion == None:  eion = np.zeros((ndens,ntemps))
    if  eele == None:  eele = np.zeros((ndens,ntemps))
    if cvion == None: cvion = np.zeros((ndens,ntemps))
    if cvele == None: cvele = np.zeros((ndens,ntemps))
    if deidn == None: deidn = np.zeros((ndens,ntemps))
    if deedn == None: deedn = np.zeros((ndens,ntemps))

    if ngroups       == None: ngroups = 1
    if opac_bounds   == None: opac_bounds = (0.0,1.0)
    if rosseland     == None: rosseland = np.zeros((ndens,ntemps,ngroups))
    if planck_absorb == None: planck_absorb = np.zeros((ndens,ntemps,ngroups))
    if planck_emiss  == None: planck_emiss = np.zeros((ndens,ntemps,ngroups))

    # Write the header:
    f = open(fn,'w')
    f.write("%10i%10i\n" % (ntemps,ndens))
    f.write(" atomic #s of gases: ")
    for z in zvals: f.write("%10i" % z)
    f.write("\n relative fractions: ")
    for frac in fracs: f.write("%10.2E" % frac)
    f.write("\n")

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
            for jt in xrange(ntemps):
                count += 1

                f.write("%s" % convert(var[jd,jt,g]))
            
                if count == 4:
                    count = 0
                    f.write("\n")

    if count != 0: f.write("\n")
        
    f.write("%12i\n" % ngroups)

    write_block(temps)
    write_block(numDens)

    write_block(zbar.flatten())
    
    write_block(dzdt.flatten())
    write_block(pion.flatten()*ERG_TO_JOULE)
    write_block(pele.flatten()*ERG_TO_JOULE)
    write_block(dpidt.flatten()*ERG_TO_JOULE)
    write_block(dpedt.flatten()*ERG_TO_JOULE)
    write_block(eion.flatten()*ERG_TO_JOULE)
    write_block(eele.flatten()*ERG_TO_JOULE)
    write_block(cvion.flatten()*ERG_TO_JOULE)
    write_block(cvele.flatten()*ERG_TO_JOULE)
    write_block(deidn.flatten()*ERG_TO_JOULE)
    write_block(deedn.flatten()*ERG_TO_JOULE)

    write_block(opac_bounds)
    write_opac_block(rosseland)
    write_opac_block(planck_absorb)
    write_opac_block(planck_emiss)


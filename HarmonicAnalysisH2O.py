import numpy as np
import numpy.linalg as la
import scipy.linalg as sla
import subprocess as sub
import math as math
import itertools as itt
import matplotlib.pyplot as plt

class Constants:
    """Helps me keep track of constants and conversions,
     originally written by Mark."""
    atomic_units = {
        "wavenumbers" : 4.55634e-6,
        "angstroms" : 1/0.529177,
        "amu" : 1.000000000000000000/6.02213670000e23/9.10938970000e-28   #1822.88839  g/mol -> a.u.
    }

    masses = {
        "H" : ( 1.00782503223, "amu"),
        "O" : (15.99491561957, "amu"),
        "D" : (2.0141017778,"amu")
    }
    @classmethod
    def convert(cls, val, unit, to_AU = True):
        vv = cls.atomic_units[unit]
        return (val * vv) if to_AU else (val / vv)

    @classmethod
    def mass(cls, atom, to_AU = True):
        m = cls.masses[atom]
        if to_AU:
            m = cls.convert(*m)
        return m


def partridgePot(cds):
    """Calls executable calc_h2o_pot, which takes in the file hoh_coord.dat and saves the energies to hoh_pot.dat
    hoh_coord.dat has the form
    nGeoms
    hx hy hz
    hx hy hz
    ox oy oz
    h'x h'y h'z
    h'x h'y h'z
    o'x o'y o'z
    h''x h''y h''z
    h''x h''y h''z
    o''x o''y o''z
    ...
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Where nGeoms is an integer corresponding to the number of geometries you pass in
    """
    np.savetxt("PES0/hoh_coord.dat", cds.reshape(cds.shape[0] * cds.shape[1], cds.shape[2]), header=str(len(cds)),
               comments="")
    sub.run('./calc_h2o_pot', cwd='PES0')
    return np.loadtxt("PES0/hoh_pot.dat")-(  -1.9109019308531233E-006)

def genStencil(cds,dispTup,dx,dim):
    if dim == 1:
        """Generates 5-point 1D stencil for finite difference"""
        atm = dispTup[0] #atom of interest
        cd = dispTup[1] #x,y, or z
        stnShape = np.concatenate(([5], np.shape(cds)))
        stencilCds = np.zeros(stnShape)
        dxOrder = [-2*dx, -dx, 0, dx,2*dx]
        ct = 0
        for disp in dxOrder:
            dxAr = np.zeros(cds.shape)
            dxAr[atm, cd] += disp
            stencilCds[ct] = cds + dxAr
            ct+=1
        vs = partridgePot(stencilCds)
        return vs
    elif dim  == 2:
        """Generates 9-point 2D stencil for finite difference"""
        atm1 = dispTup[0]
        cd1 = dispTup[1]
        atm2 = dispTup[2]
        cd2 = dispTup[3]
        stnShape = np.concatenate(([9], np.shape(cds)))
        stencilCds = np.zeros(stnShape)
        dxOrder = [-dx, 0, dx]
        ct = 0
        for disp in dxOrder:
            for disp2 in dxOrder:
                dxAr = np.zeros(cds.shape)
                dxAr[atm1, cd1] += disp
                dxAr[atm2, cd2] += disp2
                stencilCds[ct] = cds + dxAr
                ct += 1
        vs = partridgePot(stencilCds)
        return vs.reshape(3, 3)


def finiteDiff(stencil,dx,dim):
    """Weights are well known, but taken from http://www.mathematik.uni-dortmund.de/~kuzmin/cfdintro/lecture4.pdf"""
    if dim == 1: #second derivative
        #5 point stencil
        wts = np.array([
            [-1.0, 16.0, -30.0, 16.0, -1.0],
        ]) / (12.0*dx**2)
        der = stencil.dot(wts[0])
        return der
    elif dim == 2:
        # 9 point stencil
        wts = np.array([
            [1., 0., -1.],
            [1., 0., -1.]
        ]) / (dx * 2.)
        der = (stencil.dot(wts[0])).dot(wts[1])
        return der

def genHess(cds,nEls,nAtoms,dx):
    hess = np.zeros((nEls,nEls))

    ###### Off Diagonals #############
    # tupz = [[ (math.floor(x / nAtoms), x % nAtoms) for x in y ] for y in itt.combinations(range(nAtoms*nAtoms), 2)] #indexing nonsense
    # atmInfo = np.array(tupz).reshape(len(tupz),4) #indexing nonsense
    d=np.ravel_multi_index(np.array((np.triu_indices(nAtoms*3, 1))),[nAtoms*3,nAtoms*3]) #indexing nonsense
    atmInfo = np.array((np.unravel_index(d, (nAtoms, 3, nAtoms, 3)))).T #indexing nonsense (getting the upper indices of a 4D array for all my displacements)

    #atmInfo has information (cd1,atm1,cd2,atm2) which tells the stencil function which atom and dimension to displace

    offDiags = list(itt.combinations(np.arange(nEls),2)) #indices of upper triangle of hessian
    for k in range(len(offDiags)):
        stencil = genStencil(cds,atmInfo[k],dx,dim=2)
        hess[offDiags[k]] = finiteDiff(stencil,dx,dim=2)
    hess = hess+hess.T #hessian is symmetric matrix

    ###### On Diagonals #############
    atmInfo = list(itt.product(np.arange(nAtoms), repeat=2))
    for onDiags in range(nEls):
        stencil1D = genStencil(cds,atmInfo[onDiags],dx,dim=1)
        hess[onDiags,onDiags]=finiteDiff(stencil1D,dx,dim=1)
    return hess

if __name__ == '__main__':
    dx = 1.e-3
    """Everything is in Atomic Units going into generating the Hessian."""
    geom = Constants.convert(np.array(
        [[0.9578400,0.0000000,0.0000000],
        [-0.2399535,0.9272970,0.0000000],
        [0.0000000,0.0000000,0.0000000]]
    ),"angstroms",to_AU=True) #To Bohr from angstroms
    atoms = ["H","H","O"]
    nEls = 3*len(atoms)
    hessian = genHess(geom,nEls,len(atoms),dx)
    waterMasses = np.array([ Constants.mass(a) for a in atoms ])
    WaterMassesDup = np.repeat(waterMasses,3)

    #Two flavors of Diagonalizing the Hessian. Both are correct, comment one out if you like:

    #Calculate G Matrix
    GmatI = np.diag(WaterMassesDup) #G-matrix is 1/mass on diagonal for cartesian coordinates, so doing this way I am doing (1/(1/G))
    freqs,normalModes = sla.eigh(hessian,GmatI)
    freqsCM = Constants.convert(np.sqrt(freqs),'wavenumbers',to_AU=False) #To wavenumbers

    #Manually mass weight - how people think about it (divide rows and columns by inverse square root of mass of diagonal element)
    hessianMW = hessian / np.sqrt(WaterMassesDup) / np.sqrt(WaterMassesDup[:, np.newaxis])
    freqs_2, normalModes2 = la.eigh(hessianMW)
    freqsCM_2 = Constants.convert(np.sqrt(freqs_2), 'wavenumbers', to_AU=False)
    print(freqsCM) #should have 3 negative frequencies. and 3 others close to zero.
    print(freqsCM_2) #should be the same as freqsCM, should have 3 negative frequencies. and 3 others close to zero.
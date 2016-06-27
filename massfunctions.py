# massfunctions.py: A collection of functions which convert between mass and
#                   cumulative number density
# Note: all masses and number densities assumed to be logarithmic
# Sarah Wellons 6/2016+

import numpy as np
import torrey_cmf
import scipy.interpolate as interp
from scipy.optimize import newton

# ------- Illustris mass functions -------- #
def getnum_illustris(M, z):
    """
    Converts stellar mass to number density at the given redshift
    using the Illustris mass functions.

    Parameters
    ==========
    M : Stellar mass in units of log(Msun).  May be a single value or an array.
    z : Redshift

    Returns
    =======
    N : Comoving cumulative number density in units of log(Mpc^-3), same dimensions as M.

    """
    tc = torrey_cmf.number_density()
    return tc.cmf_fit(M,z)

def getmass_illustris(N, z):
    """
    Converts number density to stellar mass at the given redshift
    using the Illustris mass functions.

    Parameters
    ==========
    N : Comoving cumulative number density in units of log(Mpc^-3).  May be a single value or an array.
    z : Redshift

    Returns
    =======
    mass : Stellar mass in units of log(Msun), same dimensions as N.

    """
    tc = torrey_cmf.number_density()
    
    if isinstance(N, np.ndarray) or isinstance(N, list):
        mass = np.zeros([len(N)])
        for i, elem in enumerate(N):
            mass[i] = tc.mass_from_density(elem, z)
    else:
        mass = tc.mass_from_density(N, z)

    return mass

# ------- ZFOURGE mass functions ------- #
# Fits from Tomczak et al 2014 (ApJ 783:85)
# Note: If given a redshift outside the fit range (e.g. z < 0.35 or z > 2.75), these functions
#       will cheerfully give you the values on the nearest boundary.  Caveat emptor.
# Note2: These functions require the mpmath module.
def getnum_zfourge(M, z, target=0):
    """
    Converts stellar mass to number density at the given redshift
    using the ZFOURGE mass functions.

    Parameters
    ==========
    M : Stellar mass in units of log(Msun).  May be a single value or an array.
    z : Redshift
    target : (Optional) Used by getmass_zfourge to invert the mass function.

    Returns
    =======
    N : Comoving cumulative number density in units of log(Mpc^-3), same dimensions as M.
    """
    if isinstance(M, np.ndarray) or isinstance(M, list):
        N = np.zeros([len(M)])
        for i, elem in enumerate(M):
            N[i] = getnum_zfourge_single(elem,z)
    else:
        N = getnum_zfourge_single(M,z)

    return N - target

def getnum_zfourge_single(M, z):
    """
    Converts stellar mass to number density at the given redshift using the ZFOURGE mass 
    functions by linearly interpolating the integrated fits between neighboring redshift bins.

    Parameters
    ==========
    M : Stellar mass in units of log(Msun), must be a single value (not an array).
    z : Redshift

    Returns
    =======
    N : Comoving cumulative number density in units of log(Mpc^-3)
    """
    from mpmath import gammainc

    par1, par2 = zfourgeparams(z)
    x = 10.**(M-par1[1])
    g1 = gammainc(par1[2]+1,a=x)
    g2 = gammainc(par1[4]+1,a=x)
    N1 = np.log10(10.**(par1[3])*float(g1) + 10.**(par1[5])*float(g2))

    x = 10.**(M-par2[1])
    g1 = gammainc(par2[2]+1,a=x)
    g2 = gammainc(par2[4]+1,a=x)
    N2 = np.log10(10.**(par2[3])*float(g1) + 10.**(par2[5])*float(g2))

    return (N1*(par2[0]-z)+N2*(z-par1[0]))/(par2[0]-par1[0])

def getmass_zfourge(N, z):
    """
    Converts number density to stellar mass at the given redshift
    using the ZFOURGE mass functions.

    Parameters
    ==========
    N : Comoving cumulative number density in units of log(Mpc^-3).  May be a single value or an array.
    z : Redshift

    Returns
    =======
    mass : Stellar mass in units of log(Msun), same dimensions as N.
    """
    if isinstance(N, np.ndarray) or isinstance(N, list):
        mass = np.zeros([len(N)])
        for i, elem in enumerate(N):
            mass[i] = newton(getnum_zfourge, 10., args=(z,elem))
    else:
        mass = newton(getnum_zfourge, 10., args=(z,N))

    return mass

# Calculates the quenched fraction for the given galaxy mass and redshift
# by comparing the star-forming and quiescent mass functions
def quenchedfrac_zfourge(M, z):
    par_q1, par_q2 = zfourgeparams(z, type='quiescent')
    x_q1 = 10.**(M-par_q1[1])
    dn_q1 = np.log10(np.log(10)*np.exp(-1.*x_q1)*x_q1*(10.**par_q1[3]*x_q1**(par_q1[2]) + 10.**(par_q1[5])*x_q1**(par_q1[4])))
    x_q2 = 10.**(M-par_q2[1])
    dn_q2 = np.log10(np.log(10)*np.exp(-1.*x_q2)*x_q2*(10.**par_q2[3]*x_q2**(par_q2[2]) + 10.**(par_q2[5])*x_q2**(par_q2[4])))

    par_sf1, par_sf2 = zfourgeparams(z, type='star-forming')
    x_sf1 = 10.**(M-par_sf1[1])
    dn_sf1 = np.log10(np.log(10)*np.exp(-1.*x_sf1)*x_sf1*(10.**par_sf1[3]*x_sf1**(par_sf1[2]) + 10.**(par_sf1[5])*x_sf1**(par_sf1[4])))
    x_sf2 = 10.**(M-par_sf2[1])
    dn_sf2 = np.log10(np.log(10)*np.exp(-1.*x_sf2)*x_sf2*(10.**par_sf2[3]*x_sf2**(par_sf2[2]) + 10.**(par_sf2[5])*x_sf2**(par_sf2[4])))

    fq1 = 10.**dn_q1/(10.**dn_q1+10.**dn_sf1)
    fq2 = 10.**dn_q2/(10.**dn_q2+10.**dn_sf2)

    return (fq1*(par_q2[0]-z)+fq2*(z-par_q1[0]))/(par_q2[0]-par_q1[0])

# Returns the double-Schechter fit parameters z, log(M*), alpha1, log(P1), alpha2, and log(P2)
# of the neighboring redshift bins
def zfourgeparams(z, type='total'):
    zarr = np.array([0.35, 0.625, 0.875, 1.125, 1.375, 1.75, 2.25, 2.75]) 
    if type == 'star-forming':
        Mchar = np.array([10.59, 10.65, 10.56, 10.44, 10.69, 10.59, 10.58, 10.61])
        a1 = np.array([-1.08, -0.97, -0.46, 0.53, -0.55, 0.75, 2.06, 2.36])
        P1 = np.array([-2.67, -2.97, -2.81, -2.98, -3.04, -3.37, -4.30, -4.95])
        a2 = np.array([-2.00, -1.58, -1.61, -1.44, -1.62, -1.47, -1.38, -1.67])
        P2 = np.array([-4.46, -3.34, -3.36, -3.11, -3.59, -3.28, -3.28, -3.71])
    elif type == 'quiescent':
        Mchar = np.array([10.75, 10.68, 10.63, 10.63, 10.49, 10.77, 10.69, 9.95]) 
        a1 = np.array([-0.47, -0.10, 0.04, 0.11, 0.85, -0.19, -0.37, -0.62]) 
        P1 = np.array([-2.76, -2.67, -2.81, -3.03, -3.36, -3.41, -3.59, -4.22])
        a2 = np.array([-1.97, -1.69, -1.51, -1.57, -0.54, -0.18, -3.07, 2.51])
        P2 = np.array([-5.21, -4.29, -4.40, -4.80, -3.72, -3.91, -6.95, -4.51])
    else:   # Total
        Mchar = np.array([10.78, 10.70, 10.66, 10.54, 10.61, 10.74, 10.69, 10.74])
        a1 = np.array([-0.98, -0.39, -0.37, 0.30 , -0.12, 0.04 , 1.03 , 1.62])
        P1 = np.array([-2.54, -2.55, -2.56, -2.72, -2.78, -3.05, -3.80, -4.54])
        a2 = np.array([-1.90, -1.53, -1.61, -1.45, -1.56, -1.49, -1.33, -1.57])
        P2 = np.array([-4.29, -3.15, -3.39, -3.17, -3.43, -3.38, -3.26, -3.69])

    if z < zarr[0]:
        return [0, Mchar[0], a1[0], P1[0], a2[0], P2[0]], [zarr[0], Mchar[0], a1[0], P1[0], a2[0], P2[0]]
    elif z > zarr[-1]:
        return [zarr[-1], Mchar[-1], a1[-1], P1[-1], a2[-1], P2[-1]], [50., Mchar[-1], a1[-1], P1[-1], a2[-1], P2[-1]]
    else:
        i = np.argmax(zarr > z) - 1

    return [zarr[i], Mchar[i], a1[i], P1[i], a2[i], P2[i]], [zarr[i+1], Mchar[i+1], a1[i+1], P1[i+1], a2[i+1], P2[i+1]]

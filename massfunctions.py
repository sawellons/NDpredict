# massfunctions.py: A collection of functions which convert between mass and
#                   cumulative number density
# Note: all masses and number densities assumed to be logarithmic
# Several functions require the mpmath package.
# Sarah Wellons 6/2016+

import numpy as np
import torrey_cmf
import scipy.interpolate as interp
from scipy.optimize import newton

def getnum(M, z, massfunc='zfourge', target=0):
    """
    Converts stellar mass to number density at the given redshift using the given 
    mass function.  Note: No checks are performed to ensure that the mass function 
    is well-defined at the given parameters; it is incumbent upon the user to make 
    an appropriate choice of mass function.

    Parameters
    ==========
    M : Stellar mass in units of log(Msun).  May be a single value or an array.
    z : Redshift
    massfunc : Keyword for desired mass function.  Currently available keywords 
         include ['illustris', 'zfourge', 'muzzin', 'ilbert', 'liwhite'].
    target : Ignore - Internal parameter used for inverting mass function.

    Returns
    =======
    N : Comoving cumulative number density in units of log(Mpc^-3), same dimensions as M.

    """

    if massfunc == 'illustris':
        tc = torrey_cmf.number_density()
        return tc.cmf_fit(M,z)
    else:
        if massfunc == 'zfourge': mf = getnum_zfourge
        elif massfunc == 'muzzin': mf = getnum_muzzin
        elif massfunc == 'ilbert': mf = getnum_ilbert
        elif massfunc == 'liwhite': mf = getnum_liwhite
        else: raise ValueError("Unrecognized mass function.  Available keywords include ['illustris', 'zfourge', 'muzzin', 'ilbert', 'liwhite'].")

    if isinstance(M, np.ndarray) or isinstance(M, list):
        N = np.zeros([len(M)])
        for i, elem in enumerate(M):
            N[i] = mf(elem,z)
    else:
        N = mf(M,z)

    return N - target

def getmass(N, z, massfunc='zfourge'):
    """
    Converts number density to stellar mass at the given redshift using the the given 
    mass function.  Note: No checks are performed to ensure that the mass function 
    is well-defined at the given parameters; it is incumbent upon the user to make 
    an appropriate choice of mass function.

    Parameters
    ==========
    N : Comoving cumulative number density in units of log(Mpc^-3).  May be a single value or an array.
    z : Redshift
    massfunc : Keyword for desired mass function.  Currently available keywords 
         include ['illustris', 'zfourge', 'muzzin', 'ilbert', 'liwhite'].

    Returns
    =======
    mass : Stellar mass in units of log(Msun), same dimensions as N.

    """

    if massfunc == 'illustris':
        tc = torrey_cmf.number_density()

        if isinstance(N, np.ndarray) or isinstance(N, list):
            mass = np.zeros([len(N)])
            for i, elem in enumerate(N):
                mass[i] = tc.mass_from_density(elem, z)
        else:
            mass = tc.mass_from_density(N, z)

        return mass
    else:
        if massfunc == 'zfourge': mf = getnum_zfourge
        elif massfunc == 'muzzin': mf = getnum_muzzin
        elif massfunc == 'ilbert': mf = getnum_ilbert
        elif massfunc == 'liwhite': mf = getnum_liwhite
        else: raise ValueError("Unrecognized mass function.  Available keywords include ['illustris', 'zfourge', 'muzzin', 'ilbert', 'liwhite'].")

    if isinstance(N, np.ndarray) or isinstance(N, list):
        mass = np.zeros([len(N)])
        for i, elem in enumerate(N):
            mass[i] = newton(getnum, 10., args=(z,massfunc,elem))
    else:
        mass = newton(getnum, 10., args=(z,massfunc,N))

    return mass


# ------- COSMOS/Ultravista mass functions ------- #
# Fits from Ilbert et al. 2013 (A&A 556:55)
def getnum_ilbert(M, z):
    from mpmath import gammainc

    par1, par2 = ilbertparams(z)
    x = 10.**(M-par1[1])
    g1 = gammainc(par1[2]+1,a=x)
    g2 = gammainc(par1[4]+1,a=x)
    N1 = np.log10(par1[3]*float(g1) + par1[5]*float(g2))

    x = 10.**(M-par2[1])
    g1 = gammainc(par2[2]+1,a=x)
    g2 = gammainc(par2[4]+1,a=x)
    N2 = np.log10(par2[3]*float(g1) + par2[5]*float(g2))

    return (N1*(par2[0]-z)+N2*(z-par1[0]))/(par2[0]-par1[0])

def ilbertparams(z, type='total'):
    zarr = np.array([0.35, 0.65, 0.95, 1.3, 1.75, 2.25, 2.75, 3.5])
    Mchar = np.array([10.88, 11.03, 10.87, 10.71, 10.74, 10.74, 10.76, 10.74])
    P1 = np.array([1.68, 1.22, 2.03, 1.35, 0.88, 0.62, 0.26, 0.03])*1.e-3
    a1 = np.array([-0.69, -1., -0.52, -0.08, -0.24, -0.22, -0.15, 0.95])
    P2 = np.array([0.77, 0.16, 0.29, 0.67, 0.33, 0.15, 0.14, 0.09])*1.e-3
    a2 = np.array([-1.42, -1.64, -1.62, -1.46, -1.6, -1.6, -1.6, -1.6])

    if z < zarr[0]:
        return [0, Mchar[0], a1[0], P1[0], a2[0], P2[0]], [zarr[0], Mchar[0], a1[0], P1[0], a2[0], P2[0]]
    elif z > zarr[-1]:
        return [zarr[-1], Mchar[-1], a1[-1], P1[-1], a2[-1], P2[-1]], [50., Mchar[-1], a1[-1], P1[-1], a2[-1], P2[-1]]
    else:
        i = np.argmax(zarr > z) - 1

    return [zarr[i], Mchar[i], a1[i], P1[i], a2[i], P2[i]], [zarr[i+1], Mchar[i+1], a1[i+1], P1[i+1], a2[i+1], P2[i+1]]

# Fits from Muzzin et al. 2013 (ApJ 777:18)
def getnum_muzzin(M, z):
    from mpmath import gammainc

    par1, par2 = muzzinparams(z)
    x = 10.**(M-par1[1])
    g1 = gammainc(par1[2]+1,a=x)
    g2 = gammainc(par1[4]+1,a=x)
    N1 = np.log10(par1[3]*float(g1) + par1[5]*float(g2))

    x = 10.**(M-par2[1])
    g1 = gammainc(par2[2]+1,a=x)
    g2 = gammainc(par2[4]+1,a=x)
    N2 = np.log10(par2[3]*float(g1) + par2[5]*float(g2))

    return (N1*(par2[0]-z)+N2*(z-par1[0]))/(par2[0]-par1[0])

def muzzinparams(z, type='total'):
    zarr = np.array([0.35, 0.75, 1.25, 1.75, 2.25, 2.75, 3.5]) 
    # Alpha free
    # Mchar = np.array([10.97, 11., 10.87, 10.81, 10.81, 11.03, 11.49])
    # P1 = np.array([16.27, 16.25, 13.91, 10.13, 4.79, 1.93, 0.09])*1.e-4
    # a1 = np.array([-0.53, -1.17, -1.02, -0.86, -0.55, -1.01, -1.45])
    # P2 = np.array([9.47, 0., 0., 0., 0., 0., 0.])*1.e-4
    # a2 = np.array([-1.37, -1.2, -1.2, -1.2, -1.2, -1.2, -1.2])
    # Alpha fixed
    Mchar = np.array([10.97, 11.04, 10.99, 10.96, 11., 11.09, 11.4]) 
    P1 = np.array([16.27, 14.48, 9.30, 6.33, 2.94, 1.66, 0.13])*1.e-4
    a1 = np.array([-0.53, -1.2, -1.2, -1.2, -1.2, -1.2, -1.2])
    P2 = np.array([9.47, 0., 0., 0., 0., 0., 0.])*1.e-4
    a2 = np.array([-1.37, -1.2, -1.2, -1.2, -1.2, -1.2, -1.2])

    if z < zarr[0]:
        return [0, Mchar[0], a1[0], P1[0], a2[0], P2[0]], [zarr[0], Mchar[0], a1[0], P1[0], a2[0], P2[0]]
    elif z > zarr[-1]:
        return [zarr[-1], Mchar[-1], a1[-1], P1[-1], a2[-1], P2[-1]], [50., Mchar[-1], a1[-1], P1[-1], a2[-1], P2[-1]]
    else:
        i = np.argmax(zarr > z) - 1

    return [zarr[i], Mchar[i], a1[i], P1[i], a2[i], P2[i]], [zarr[i+1], Mchar[i+1], a1[i+1], P1[i+1], a2[i+1], P2[i+1]]


# ------- Low-z mass function from Li & White 2009 (MNRAS 398:2977) ------ #
# Note: Only applies to a single redshift (z ~ 0.1)!
def getnum_liwhite(M, z):
    from mpmath import gammainc

    h = 0.7
    # log(M*), alpha, Phi*,
    par_high = [10.71 - np.log10(h**2), -1.99, 0.0044*h**3] 
    par_mid = [10.37 - np.log10(h**2), -0.9, 0.0132*h**3] 
    par_low = [9.61 - np.log10(h**2), -1.13, 0.0146*h**3] 

    # Pre-tabulated integrals
    hightot = 0.000280917465551   # Total contribution from high-mass piece
    subtractmid = 0.000245951678324549    # Value of middle Schechter fn integrated down to first break
    midtot = 0.00748165061711    # Total contribution from mid-mass piece
    subtractlow = 0.00266964412065

    if M > 10.67 - np.log10(h**2):
        x = 10.**(M-par_high[0])
        g = gammainc(par_high[1]+1,a=x)
        return np.log10(par_high[2]*float(g))
    elif M > 9.33 - np.log10(h**2):
        x = 10.**(M-par_mid[0])
        g = gammainc(par_mid[1]+1,a=x)
        return np.log10(par_mid[2]*float(g)-subtractmid+hightot)
    else:
        x = 10.**(M-par_low[0])
        g = gammainc(par_low[1]+1,a=x)
        return np.log10(par_low[2]*float(g)-subtractlow+midtot+hightot)

# ------- ZFOURGE mass functions ------- #
# Fits from Tomczak et al 2014 (ApJ 783:85)
def getnum_zfourge(M, z):
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

# Returns the double-Schechter fit parameters z, log(M*), alpha1, log(P1), alpha2, and log(P2)
# of the neighboring redshift bins
def zfourgeparams(z, type='total'):
    zarr = np.array([0.35, 0.625, 0.875, 1.125, 1.375, 1.75, 2.25, 2.75]) 
    # Double-Schechter at z > 2
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

    # # Single-Schechter at z > 2
    # if type == 'star-forming':
    #     Mchar = np.array([10.59, 10.65, 10.56, 10.44, 10.69, 10.59, 11.28, 11.49])
    #     a1 = np.array([-1.08, -0.97, -0.46, 0.53, -0.55, 0.75, -1.6, -1.93])
    #     P1 = np.array([-2.67, -2.97, -2.81, -2.98, -3.04, -3.37, -3.96, -4.82])
    #     a2 = np.array([-2.00, -1.58, -1.61, -1.44, -1.62, -1.47, -1.38, -1.67])
    #     P2 = np.array([-4.46, -3.34, -3.36, -3.11, -3.59, -3.28, -50., -50.])
    # elif type == 'quiescent':
    #     Mchar = np.array([10.75, 10.68, 10.63, 10.63, 10.49, 10.77, 10.73, 10.65]) 
    #     a1 = np.array([-0.47, -0.10, 0.04, 0.11, 0.85, -0.19, -0.49, -0.43]) 
    #     P1 = np.array([-2.76, -2.67, -2.81, -3.03, -3.36, -3.41, -3.63, -3.92])
    #     a2 = np.array([-1.97, -1.69, -1.51, -1.57, -0.54, -0.18, -3.07, 2.51])
    #     P2 = np.array([-5.21, -4.29, -4.40, -4.80, -3.72, -3.91, -50., -50.])
    # else:   # Total
    #     Mchar = np.array([10.78, 10.70, 10.66, 10.54, 10.61, 10.74, 11.13, 11.35])
    #     a1 = np.array([-0.98, -0.39, -0.37, 0.30 , -0.12, 0.04, -1.43, -1.74])
    #     P1 = np.array([-2.54, -2.55, -2.56, -2.72, -2.78, -3.05, -3.59, -4.36])
    #     a2 = np.array([-1.90, -1.53, -1.61, -1.45, -1.56, -1.49, -1.33, -1.57])
    #     P2 = np.array([-4.29, -3.15, -3.39, -3.17, -3.43, -3.38, -50., -50.])

    if z < zarr[0]:
        return [0, Mchar[0], a1[0], P1[0], a2[0], P2[0]], [zarr[0], Mchar[0], a1[0], P1[0], a2[0], P2[0]]
    elif z > zarr[-1]:
        return [zarr[-1], Mchar[-1], a1[-1], P1[-1], a2[-1], P2[-1]], [50., Mchar[-1], a1[-1], P1[-1], a2[-1], P2[-1]]
    else:
        i = np.argmax(zarr > z) - 1

    return [zarr[i], Mchar[i], a1[i], P1[i], a2[i], P2[i]], [zarr[i+1], Mchar[i+1], a1[i+1], P1[i+1], a2[i+1], P2[i+1]]

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

# ------ OBSOLETE, left in for backwards-compatibility ------ #
def getnum_illustris(M, z):
    tc = torrey_cmf.number_density()
    return tc.cmf_fit(M,z)

def getmass_illustris(N, z):
    tc = torrey_cmf.number_density()
    
    if isinstance(N, np.ndarray) or isinstance(N, list):
        mass = np.zeros([len(N)])
        for i, elem in enumerate(N):
            mass[i] = tc.mass_from_density(elem, z)
    else:
        mass = tc.mass_from_density(N, z)

    return mass

def getmass_zfourge(N, z):
    if isinstance(N, np.ndarray) or isinstance(N, list):
        mass = np.zeros([len(N)])
        for i, elem in enumerate(N):
            mass[i] = newton(getnum_zfourge, 10., args=(z,elem))
    else:
        mass = newton(getnum_zfourge, 10., args=(z,N))

    return mass


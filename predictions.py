# predictions.py: Predict progenitors'/descendants' stellar masses and assign probabilities
# Note: all number densities, stellar masses assumed to be logarithmic
# Sarah Wellons 6/2016+

import numpy as np
from .Nevol import *
from .massfunctions import *

def assign_probabilities(M0, z0, zf, M_sample, volume, Nfunc=evolvingN, sigmafunc=sigmaN, massfunc='zfourge'):
    """
    For every galaxy in a sample at redshift zf, find the probability that it is the progenitor/descendant 
    of a galaxy with mass M0 at redshift z0.
    [Note: This is NOT equivalent to the probability that the mass of the descendant/progenitor of the
    the galaxy is M0!  P(Mf|M0) != P(M0|Mf), see Torrey et al 2016.]

    Parameters
    ==========
    M0 : Initial stellar mass in units of log(Msun)
    z0 : Initial redshift
    zf : Final redshift, at which the prediction will be made
    M_sample : The masses of the galaxies in the sample at zf in units of log(Msun)
    volume : The volume of the sample at zf in units of Mpc^3.  This may either be a constant value 
             or an array with the same dimensions as M_sample to allow for an effective volume 
             which is a function of stellar mass.
    Nfunc : (Optional) Function which evolves number density to another redshift.
            Must take arguments (N0, z0, zf), defaults to evolvingN.
    sigmafunc : (Optional) Function which predicts the width of the number density distribution at another redshift.
                Must take arguments (N0, z0, zf), defaults to sigmaN.
    massfunc : (Optional) Keyword for mass function used to convert from mass to number density. 
           Currently available keywords include ['illustris', 'zfourge', 'muzzin', 'ilbert', 'liwhite'].
           Defaults to 'zfourge'.

    Returns
    =======
    Array of progenitor probabilities, the same dimension as M_sample
    """
    N0 = getnum(M0, z0, massfunc=massfunc)
    N_sample = getmass(M_sample, zf, massfunc=massfunc)
    dpdlogN = newnum_distrib(N0, z0, zf, N_sample, Nfunc=Nfunc, sigmafunc=sigmafunc)
    pgal = dpdlogN / 10.**N_sample / volume
    return pgal 

def newnum_distrib(N0, z0, zf, Narr, Nfunc=evolvingN, sigmafunc=sigmaN):
    """
    Calculates the probability density dp/dlogN that the progenitor/descendant of a 
    galaxy with number density N0 at redshift z0 has number density Narr at redshift zf.

    Parameters
    ==========
    N0 : Initial number density at z0 in units of log(Mpc^-3)
    z0 : Initial redshift
    zf : Final redshift, at which the prediction will be made
    Narr : Number density at zf (in units of log(Mpc^-3)) at which to evaluate the pdf.  
           May be single value or array.
    Nfunc : (Optional) Function which evolves number density to another redshift.
            Must take arguments (N0, z0, zf), defaults to evolvingN.
    sigmafunc : (Optional) Function which predicts the width of the number density distribution at another redshift.
                Must take arguments (N0, z0, zf), defaults to sigmaN.

    Returns
    =======
    Probability density dp/dlogN corresponding to Narr at zf.
    """
    pi = 3.14159265359
    Nf = Nfunc(N0, z0, zf)
    sigma = sigmafunc(N0, z0, zf)
    return np.exp(-(Narr-Nf)**2/(2*sigma**2))/(sigma*np.sqrt(2.*pi))

def newmass(M0, z0, zf, Nfunc=evolvingN, massfunc='zfourge'):
    """
    For a galaxy population with initial stellar mass M0 at redshift z0, predicts the median 
    stellar mass at another redshift zf.

    Parameters
    ==========
    M0 : Initial stellar mass (at z0) in units of log(Msun)
    z0 : Initial redshift.  If zf > z0, z0 must equal 0.
    zf : Final redshift for which the prediction is to be made.
         May be greater than or less than z0.
    Nfunc : (Optional) Function which evolves number density to another redshift.
            Must take arguments (N0, z0, zf), defaults to evolvingN.
    massfunc : (Optional) Keyword for mass function used to convert between mass and number density. 
           Currently available keywords include ['illustris', 'zfourge', 'muzzin', 'ilbert', 'liwhite'].
           Defaults to 'zfourge'.

    Returns
    =======
    Mf : Median stellar mass (in units of log(Msun)) of the population at zf.
    """

    N0 = getnum(M0, z0, massfunc=massfunc)
    Nf = Nfunc(N0, z0, zf)
    return getmass(Nf, zf, massfunc=massfunc)

def newmass_distrib(M0, z0, zf, Marr=None, Medges=None, Nfunc=evolvingN, sigmafunc=sigmaN, massfunc='zfourge'):
    """
    Evaluates the probability that progenitors/descendants lie in the given mass bins at zf.

    Parameters
    ==========
    M0 : Initial stellar mass (at z0), logarithmic
    z0 : Initial redshift
    zf : Final redshift, at which the prediction will be made
    Marr : (Optional) Array of stellar masses central to each bin, length n
    Medges : (Optional) Array of edges of the mass bins, length n+1
             At least one of Marr or Medges must be specified.
    Nfunc : (Optional) Function which evolves number density to another redshift.
            Must take arguments (N0, z0, zf), defaults to evolvingN.
    sigmafunc : (Optional) Function which predicts the width of the number density distribution at another redshift.
                Must take arguments (N0, z0, zf), defaults to sigmaN.
    massfunc : (Optional) Keyword for mass function used to convert between mass and number density. 
           Currently available keywords include ['illustris', 'zfourge', 'muzzin', 'ilbert', 'liwhite'].
           Defaults to 'zfourge'.

    Returns
    =======
    Array of length n representing the probability that a progenitor/descendant lies in that mass bin.
    """

    # Check for z0 = zf so there's no division by sigma = 0
    if z0 == zf:
        parr = np.zeros([len(Marr)])
        Mind = np.arange(0,len(Medges))
        Mi = np.min(Mind[Medges > M0])
        parr[Mi-1] = 1./(Medges[Mi] - Medges[Mi-1])
        return parr

    if Medges is None:
        try:
            Medges = np.zeros(len(Marr)+1)
            Medges[1:-1] = (Marr[1:]+Marr[:-1])/2.
            Medges[0] = 2.*Medges[1] - Medges[2]
            Medges[-1] = 2.*Medges[-2] - Medges[-3]
        except:
            raise ValueError("Mass array is too short to define edges.  Please supply your own.")

    if Marr is None:
        try:
            Marr = (Medges[1:]+Medges[:-1])/2.
        except:
            raise ValueError("Medges must have length of at least 2.")

    Narr = getnum(Marr, zf, massfunc=massfunc)
    Nedges = getnum(Medges, zf, massfunc=massfunc)
    dN = abs(Nedges[1:] - Nedges[:-1])
    return newnum_distrib(getnum(M0, z0, massfunc=massfunc), z0, zf, Narr, Nfunc=Nfunc, sigmafunc=sigmafunc) * dN


# Nevol.py: A collection of functions which predict the median number density 
#           of a population at redshift z given an initial number density at redshift z0
# Note: all number densities assumed to be logarithmic
# Sarah Wellons 6/2016+

import numpy as np
from subprocess import check_output
import torrey_cmf
tc = torrey_cmf.number_density()

def evolvingN(N0, z0, zf, type='IllustrisCMF'):
    """
    Predicts the median number density at zf for a galaxy population with initial
    number density N0 at redshift z0, using results drawn from Illustris or Millennium.

    Parameters
    ==========
    N0 : Initial comoving cumulative number density (in units of log(Mpc^-3)) at z0.
    z0 : Initial redshift.  If zf > z0, z0 must equal 0.
    zf : Final redshift for which the prediction is to be made.  
         May be greater than or less than z0.
    type : (Optional) Keyword which indicates whether the Illustris or Millennium tracks
           are to be employed.  Defaults to Illustris.

    Returns
    =======
    Nf : Median number density of the population at zf, in units of log(Mpc^-3).

    """
    if z0 < zf:
        if abs(z0) < 0.01:
            return tc.nd_backward_fit(zf, N0, type=type)
        else:
            raise ValueError("If predicting back in time, must begin at z=0.")
    else:
        return tc.single_nd_fit(zf, z0, N0, type=type)

def milN(N0, z0, zf):
    return evolvingN(N0, z0, zf, type='MillenniumCMF')

def constantN(N0, z0, zf):
    """ Comparable to evolvingN, but for a constant number density. """
    return N0

def sigmaN(N0, z0, zf, type='IllustrisCMF'):
    """
    Predicts the width of the logarithmic number density distribution at zf for a 
    galaxy population with initial number density N0 at redshift z0, using results 
    drawn from Illustris or Millennium.

    Parameters
    ==========
    N0 : Initial comoving cumulative number density (in units of log(Mpc^-3)) at z0.
    z0 : Initial redshift.  If zf > z0, z0 must equal 0.
    zf : Final redshift for which the prediction is to be made.  
         May be greater than or less than z0.
    type : (Optional) Keyword which indicates whether the Illustris or Millennium tracks
           are to be employed.  Defaults to Illustris.

    Returns
    =======
    sigma : The standard deviation of the population's lognormal distribution 
            in number density at zf, in units of log(Mpc^-3).
    """
    if z0 < zf:
        if abs(z0) < 0.01:
            return tc.sigma_backward_fit(zf, N0, type=type)
        else:
            raise ValueError("If predicting back in time, must begin at z=0.")
    else:
        return tc.sigma_forward_fit(zf, z0, N0, type=type)

def milsigma(N0, z0, zf):
    return sigmaN(N0, z0, zf, type='MillenniumCMF')

def behrooziN(N0, z0, zf):
    """
    Predicts the median number density and +- 1 sigma countours at zf for a galaxy 
    population with initial number density N0 at redshift z0, using results drawn 
    from Bolshoi by Behroozi et al (2013).

    Parameters
    ==========
    N0 : Initial comoving cumulative number density (in units of log(Mpc^-3)) at z0.
    z0 : Initial redshift.  If zf > z0, z0 must equal 0.
    zf : Final redshift for which the prediction is to be made.  
         May be greater than or less than z0.

    Returns
    =======
    Nmed, Nminus, Nplus : [N, N-1sigma, N+1sigma] for the population at zf, in units of log(Mpc^-3).
    """

    try:
        ret = check_output(["./nd_redshift", str(z0), str(N0), str(zf)]).split()
        return float(ret[3]), float(ret[4]), float(ret[5])
    except:
        raise ValueError("Improper use of nd_redshift function.  Make sure to have the nd_redshift executable, mf_bolshoi.bin and z_tracks.bin in this directory.")

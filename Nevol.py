# Nevol.py: A collection of functions which predict the median number density 
#           of a population at redshift z given an initial number density at redshift z0
# Note: all number densities assumed to be logarithmic
# Sarah Wellons 6/2016+

import numpy as np
import torrey_cmf
tc = torrey_cmf.number_density()

def evolvingN(N0, z0, zf):
    """
    Predicts the median number density at zf for a galaxy population with initial
    number density N0 at redshift z0.

    Parameters
    ==========
    N0 : Initial comoving cumulative number density (in units of log(Mpc^-3)) at z0.
    z0 : Initial redshift.  If zf > z0, z0 must equal 0.
    zf : Final redshift for which the prediction is to be made.  
         May be greater than or less than z0.

    Returns
    =======
    Nf : Median number density of the population at zf, in units of log(Mpc^-3).

    """
    if z0 < zf:
        if abs(z0) < 0.01:
            return tc.nd_backward_fit(zf, N0)
        else:
            raise ValueError("If predicting back in time, must begin at z=0.")
    else:
        return tc.single_nd_fit(zf, z0, N0)

def constantN(N0, z0, zf):
    """ Comparable to evolvingN, but for a constant number density. """
    return N0

def sigmaN(N0, z0, zf):
    """
    Predicts the width of the logarithmic number density distribution at zf for a 
    galaxy population with initial number density N0 at redshift z0.

    Parameters
    ==========
    N0 : Initial comoving cumulative number density (in units of log(Mpc^-3)) at z0.
    z0 : Initial redshift.  If zf > z0, z0 must equal 0.
    zf : Final redshift for which the prediction is to be made.  
         May be greater than or less than z0.

    Returns
    =======
    sigma : The standard deviation of the population's lognormal distribution 
            in number density at zf, in units of log(Mpc^-3).
    """
    if z0 < zf:
        if abs(z0) < 0.01:
            return tc.sigma_backward_fit(zf, N0)
        else:
            raise ValueError("If predicting back in time, must begin at z=0.")
    else:
        return tc.sigma_forward_fit(zf, z0, N0)



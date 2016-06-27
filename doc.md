This file contains the docstrings for all user functions in NDpredict.

## Nevol.py

###evolvingN(N0, z0, zf):
Predicts the median number density at zf for a galaxy population with initial number density N0 at redshift z0.

**Parameters**
- N0 : Initial comoving cumulative number density (in units of log(Mpc^-3)) at z0.
- z0 : Initial redshift.  If zf > z0, z0 must equal 0.
- zf : Final redshift for which the prediction is to be made.  May be greater than or less than z0.

**Returns**
- Nf : Median number density of the population at zf, in units of log(Mpc^-3).

###constantN(N0, z0, zf):
    """ Comparable to evolvingN, but for a constant number density. """
    return N0

###sigmaN(N0, z0, zf):
Predicts the width of the logarithmic number density distribution at zf for a galaxy population with initial number density N0 at redshift z0.

**Parameters**
- N0 : Initial comoving cumulative number density (in units of log(Mpc^-3)) at z0.
- z0 : Initial redshift.  If zf > z0, z0 must equal 0.
- zf : Final redshift for which the prediction is to be made.  May be greater than or less than z0.

**Returns**
- sigma : The standard deviation of the population's lognormal distribution in number density at zf, in units of log(Mpc^-3).




## predictions.py

###assign_probabilities(M0, z0, zf, M_sample, volume, Nfunc=evolvingN, sigmafunc=sigmaN, MtoN=getnum_zfourge):
For every galaxy in a sample at redshift zf, find the probability that it is the progenitor/descendant of a galaxy with mass M0 at redshift z0.
[Note: This is NOT equivalent to the probability that the mass of the descendant/progenitor of the galaxy is M0!  P(Mf|M0) != P(M0|Mf), see Torrey et al 2016.]

**Parameters**
- M0 : Initial stellar mass in units of log(Msun)
- z0 : Initial redshift
- zf : Final redshift, at which the prediction will be made
- M_sample : The masses of the galaxies in the sample at zf in units of log(Msun)
- volume : The volume of the sample at zf in units of Mpc^3.  This may either be a constant value or an array with the same dimensions as M_sample to allow for an effective volume which is a function of stellar mass.
- Nfunc : (Optional) Function which evolves number density to another redshift. Must take arguments (N0, z0, zf), defaults to evolvingN.
- sigmafunc : (Optional) Function which predicts the width of the number density distribution at another redshift.  Must take arguments (N0, z0, zf), defaults to sigmaN.
- MtoN : (Optional) Function which converts from mass to number density.  Must take arguments (N0, z0, zf), defaults to getnum_zfourge.

**Returns**
- Array of progenitor probabilities, the same dimension as M_sample

###newnum_distrib(N0, z0, zf, Narr, Nfunc=evolvingN, sigmafunc=sigmaN):
Calculates the probability density dp/dlogN that the progenitor/descendant of a galaxy with number density N0 at redshift z0 has number density Narr at redshift zf.

**Parameters**
- N0 : Initial number density at z0 in units of log(Mpc^-3)
- z0 : Initial redshift
- zf : Final redshift, at which the prediction will be made
- Narr : Number density at zf (in units of log(Mpc^-3)) at which to evaluate the pdf.  May be single value or array.
- Nfunc : (Optional) Function which evolves number density to another redshift.  Must take arguments (N0, z0, zf), defaults to evolvingN.
- sigmafunc : (Optional) Function which predicts the width of the number density distribution at another redshift.  Must take arguments (N0, z0, zf), defaults to sigmaN.

**Returns**
- Probability density dp/dlogN corresponding to Narr at zf.

###newmass(M0, z0, zf, Nfunc=evolvingN, MtoN=getnum_zfourge, NtoM=getmass_zfourge):
For a galaxy population with initial stellar mass M0 at redshift z0, predicts the median stellar mass at another redshift zf.

**Parameters**
- M0 : Initial stellar mass (at z0) in units of log(Msun)
- z0 : Initial redshift.  If zf > z0, z0 must equal 0.
- zf : Final redshift for which the prediction is to be made.  May be greater than or less than z0.
- Nfunc : (Optional) Function which evolves number density to another redshift.  Must take arguments (N0, z0, zf), defaults to evolvingN.
- MtoN : (Optional) Function which converts from mass to number density.  Must take arguments (N0, z0, zf), defaults to getnum_zfourge.
- NtoM : (Optional) Function which converts from number density to mass.  Must take arguments (N0, z0, zf), defaults to getmass_zfourge.

**Returns**
- Mf : Median stellar mass (in units of log(Msun)) of the population at zf.

###newmass_distrib(M0, z0, zf, Marr=None, Medges=None, Nfunc=evolvingN, sigmafunc=sigmaN, MtoN=getnum_zfourge):
Evaluates the probability that progenitors/descendants lie in the given mass bins at zf.

**Parameters**
- M0 : Initial stellar mass (at z0), logarithmic
- z0 : Initial redshift
- zf : Final redshift, at which the prediction will be made
- Marr : (Optional) Array of stellar masses central to each bin, length n.
- Medges : (Optional) Array of edges of the mass bins, length n+1.   At least one of Marr or Medges must be specified.
- Nfunc : (Optional) Function which evolves number density to another redshift.  Must take arguments (N0, z0, zf), defaults to evolvingN.
- sigmafunc : (Optional) Function which predicts the width of the number density distribution at another redshift.  Must take arguments (N0, z0, zf), defaults to sigmaN.
- MtoN : (Optional) Function which converts from mass to number density.  Must take arguments (N0, z0, zf), defaults to getnum_zfourge.

**Returns**
- Array of length n representing the probability that a progenitor/descendant lies in that mass bin.


## massfunctions.py
 
###getnum_illustris(M, z):
Converts stellar mass to number density at the given redshift using the Illustris mass functions.

**Parameters**
- M : Stellar mass in units of log(Msun).  May be a single value or an array.
- z : Redshift

**Returns**
- N : Comoving cumulative number density in units of log(Mpc^-3), same dimensions as M.

###getmass_illustris(N, z):
Converts number density to stellar mass at the given redshift using the Illustris mass functions.

**Parameters**
- N : Comoving cumulative number density in units of log(Mpc^-3).  May be a single value or an array.
- z : Redshift

**Returns**
- mass : Stellar mass in units of log(Msun), same dimensions as N.

###getnum_zfourge(M, z):
Converts stellar mass to number density at the given redshift using the ZFOURGE mass functions.

**Parameters**
- M : Stellar mass in units of log(Msun).  May be a single value or an array.
- z : Redshift

**Returns**
- N : Comoving cumulative number density in units of log(Mpc^-3), same dimensions as M.

###def getmass_zfourge(N, z):
Converts number density to stellar mass at the given redshift using the Illustris mass functions.

**Parameters**
- N : Comoving cumulative number density in units of log(Mpc^-3).  May be a single value or an array.
- z : Redshift

**Returns**
- mass : Stellar mass in units of log(Msun), same dimensions as N.


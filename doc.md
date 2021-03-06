This file contains the docstrings for all user functions in NDpredict.

## Nevol.py

###evolvingN(N0, z0, zf, type='IllustrisCMF'):
Predicts the median number density at zf for a galaxy population with initial number density N0 at redshift z0, using results drawn from Illustris or Millennium.

**Parameters**
- N0 : Initial comoving cumulative number density (in units of log(Mpc^-3)) at z0.
- z0 : Initial redshift.  If zf > z0, z0 must equal 0.
- zf : Final redshift for which the prediction is to be made.  May be greater than or less than z0.
- type : (Optional) Keyword which indicates whether the Illustris or Millennium tracks are to be employed.  Defaults to Illustris.

**Returns**
- Nf : Median number density of the population at zf, in units of log(Mpc^-3).

###milN(N0, z0, zf)
Calls evolvingN with type='MillenniumCMF'

###constantN(N0, z0, zf):
Comparable to evolvingN, but for a constant number density. 

###sigmaN(N0, z0, zf, type='IllustrisCMF'):
Predicts the width of the logarithmic number density distribution at zf for a galaxy population with initial number density N0 at redshift z0, using results drawn from Illustris or Millennium.

**Parameters**
- N0 : Initial comoving cumulative number density (in units of log(Mpc^-3)) at z0.
- z0 : Initial redshift.  If zf > z0, z0 must equal 0.
- zf : Final redshift for which the prediction is to be made.  May be greater than or less than z0.
- type : (Optional) Keyword which indicates whether the Illustris or Millennium tracks are to be employed.  Defaults to Illustris.

**Returns**
- sigma : The standard deviation of the population's lognormal distribution in number density at zf, in units of log(Mpc^-3).

###milsigma(N0, z0, zf)
Calls sigmaN with type='MillenniumCMF'

###behrooziN(N0, z0, zf, type='IllustrisCMF'):
Predicts the median number density and +- 1 sigma contours at zf for a galaxy population with initial number density N0 at redshift z0, using results drawn from Bolshoi by Behroozi et al (2013).  Requires the nd_redshift executable, mf_bolshoi.bin and z_tracks.bin from https://code.google.com/archive/p/nd-redshift/.

**Parameters**
- N0 : Initial comoving cumulative number density (in units of log(Mpc^-3)) at z0.
- z0 : Initial redshift.  
- zf : Final redshift for which the prediction is to be made.  May be greater than or less than z0.

**Returns**
- Nmed, Nminus, Nplus : [N, N-1sigma, N+1sigma] for the population at zf, in units of log(Mpc^-3).

## predictions.py

###assign_probabilities(M0, z0, zf, M_sample, volume, Nfunc=evolvingN, sigmafunc=sigmaN, massfunc='zfourge'):
For every galaxy in a sample at redshift zf, find the probability that it is the progenitor/descendant of a galaxy with mass M0 at redshift z0.
[Note: This is NOT equivalent to the probability that the mass of the descendant/progenitor of the galaxy is M0!  P(Mf|M0) != P(M0|Mf), see Torrey et al 2017.]

**Parameters**
- M0 : Initial stellar mass in units of log(Msun)
- z0 : Initial redshift
- zf : Final redshift, at which the prediction will be made
- M_sample : The masses of the galaxies in the sample at zf in units of log(Msun)
- volume : The volume of the sample at zf in units of Mpc^3.  This may either be a constant value or an array with the same dimensions as M_sample to allow for an effective volume which is a function of stellar mass.
- Nfunc : (Optional) Function which evolves number density to another redshift. Must take arguments (N0, z0, zf), defaults to evolvingN.
- sigmafunc : (Optional) Function which predicts the width of the number density distribution at another redshift.  Must take arguments (N0, z0, zf), defaults to sigmaN.
- massfunc : (Optional) Keyword for mass function used to convert from mass to number density.  Currently available keywords include ['illustris', 'zfourge', 'muzzin', 'ilbert', 'liwhite']. Defaults to 'zfourge'.

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

###newmass(M0, z0, zf, Nfunc=evolvingN, massfunc='zfourge'):
For a galaxy population with initial stellar mass M0 at redshift z0, predicts the median stellar mass at another redshift zf.

**Parameters**
- M0 : Initial stellar mass (at z0) in units of log(Msun)
- z0 : Initial redshift.  If zf > z0, z0 must equal 0.
- zf : Final redshift for which the prediction is to be made.  May be greater than or less than z0.
- Nfunc : (Optional) Function which evolves number density to another redshift.  Must take arguments (N0, z0, zf), defaults to evolvingN.
- massfunc : (Optional) Keyword for mass function used to convert from mass to number density.  Currently available keywords include ['illustris', 'zfourge', 'muzzin', 'ilbert', 'liwhite']. Defaults to 'zfourge'.

**Returns**
- Mf : Median stellar mass (in units of log(Msun)) of the population at zf.

###newmass_distrib(M0, z0, zf, Marr=None, Medges=None, Nfunc=evolvingN, sigmafunc=sigmaN, massfunc='zfourge'):
Evaluates the probability that progenitors/descendants lie in the given mass bins at zf.

**Parameters**
- M0 : Initial stellar mass (at z0), logarithmic
- z0 : Initial redshift
- zf : Final redshift, at which the prediction will be made
- Marr : (Optional) Array of stellar masses central to each bin, length n.
- Medges : (Optional) Array of edges of the mass bins, length n+1.   At least one of Marr or Medges must be specified.
- Nfunc : (Optional) Function which evolves number density to another redshift.  Must take arguments (N0, z0, zf), defaults to evolvingN.
- sigmafunc : (Optional) Function which predicts the width of the number density distribution at another redshift.  Must take arguments (N0, z0, zf), defaults to sigmaN.
- massfunc : (Optional) Keyword for mass function used to convert from mass to number density.  Currently available keywords include ['illustris', 'zfourge', 'muzzin', 'ilbert', 'liwhite']. Defaults to 'zfourge'.

**Returns**
- Array of length n representing the probability that a progenitor/descendant lies in that mass bin.


## massfunctions.py
 
The two general functions getnum(M,z) and getmass(N,z) can be used for any of the mass functions described below using the appropriate keyword for the optional argument 'massfunc'.  Note that these are cumulative mass functions, not differential.  Supported mass functions currently include:
- 'illustris' : From Torrey et al. (2015), stellar mass functions from the Illustris simulation from 0 < z < 3.
- 'zfourge' : From Tomczak et al. (2014), stellar mass functions from the ZFOURGE survey from 0.2 < z < 3.
- 'muzzin' : From Muzzin et al. (2013), stellar mass functions from the COSMOS/UltraVISTA survey from 0.2 < z < 4.
- 'ilbert' : From Ilbert et al. (2013), stellar mass functions from the COSMOS/UltraVISTA survey from 0.2 < z < 4.
- 'liwhite': From Li & White (2009), stellar mass function from SDSS at z < 0.1
 
###getnum(M, z, massfunc='zfourge', interpdir='N'):
Converts stellar mass to number density at the given redshift using the given mass function.  Note: No checks are performed to ensure that the mass function is well-defined at the given parameters; it is incumbent upon the user to make an appropriate choice of mass function.

**Parameters**
- M : Stellar mass in units of log(Msun).  May be a single value or an array.
- z : Redshift
- massfunc : Keyword for desired mass function, as listed above.  
- interpdir : (Optional) ['N', 'M'] Indicates whether the mass functions should be interpolated across redshift in log N or log M.  Doesn't make a significant difference in most cases.

**Returns**
- N : Comoving cumulative number density in units of log(Mpc^-3), same dimensions as M.

###getmass(N, z, massfunc='zfourge'):
Converts number density to stellar mass at the given redshift.

**Parameters**
- N : Comoving cumulative number density in units of log(Mpc^-3).  May be a single value or an array.
- z : Redshift
- massfunc : Keyword for desired mass function, as listed above.  
- interpdir : (Optional) ['N', 'M'] Indicates whether the mass functions should be interpolated across redshift in log N or log M.  Doesn't make a significant difference in most cases.

**Returns**
- mass : Stellar mass in units of log(Msun), same dimensions as N.


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

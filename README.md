# NDpredict

NDpredict may be used with observational or simulation data to predict the motion of galaxy populations through cumulative comoving number density space, and assign progenitor/descendant probabilities to a galaxy sample.

## Requirements

If you employ this package in your research, please cite 
- Wellons and Torrey 2016 (arxiv:1606.07815)
- Torrey et al 2016 (arxiv:1606.07271)

To use the number density evolution functions or the mass functions from Illustris, you will need the torrey_cmf package (https://github.com/ptorrey/torrey_cmf).

To use any of the observational mass functions, you will need the mpmath package (mpmath.org).

## Functionality

All functions are documented in detail in doc.md.

### Predicting number density / stellar mass evolution

The prediction of number densities can be done with the functions evolvingN and sigmaN in Nevol.py, e.g.
```
import NDpredict as ndp
z0 = 2.
zf = 0.
N0 = -3.   # All number densities are assumed to be logarithmic
print(ndp.evolvingN(N0, z0, zf))
>> -3.0813739999997916
print(ndp.sigmaN(N0, z0, zf))
>> 0.34528558000000015
```

To work only in terms of stellar mass, one may use the functions newmass or newmass_distrib in predictions.py.  A mass function must be specified to use these functions (and defaults to the ZFOURGE functions), e.g.
```
import NDpredict as ndp
z0 = 2.
zf = 0.
M0 = 10.5   # All stellar masses are assumed to be logarithmic
print(ndp.newmass(M0, z0, zf, massfunc='zfourge')  # Median stellar mass at zf
>> 10.845628052372005
Marray = numpy.linspace(10.5, 11.5, 10)
# Probability that a descendant falls into these mass bins
print(ndp.newmass_distrib(M0, z0, zf, Medges=Marray, massfunc='zfourge'))
>> array([  7.06694540e-02,   1.26868211e-01,   1.99301323e-01,
         2.45146482e-01,   1.97076762e-01,   7.72581613e-02,
         9.18849092e-03,   1.53544788e-04,   1.03212430e-07])
```


### Mass Functions

The available mass functions can be found in massfunctions.py and currently include:
- Illustris (Torrey et al 2015)
- ZFOURGE (Tomczak et al 2014)
- Muzzin et al 2013
- Ilbert et al 2013
- Li & White 2009

If you would like to use another set of mass functions, please email swellons@cfa.harvard.edu to have them added.

### Assigning progenitor/descendant probabilities

The assignment of progenitor/descendant probabilities to galaxies in a sample can be done with the function assign_probabilities in predictions.py.  Distributions in progenitor/descendant galaxy properties can then be predicted by using the probabilities as weights and summing over the entire galaxy population, e.g.
```
import NDpredict as ndp
sample_masses = get_sample() # Array of stellar masses of all galaxies in the sample
vol = 100.^3  # Sample volume in Mpc^3
z0 = 2.
zf = 0.
M0 = 10.5  
probs = ndp.assign_probabilities(M0, z0, zf, numpy.log10(sample_masses), vol, massfunc='zfourge')
# Mean descendant mass
Mavg = numpy.average(sample_masses, weights=probs)
```

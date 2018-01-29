# Parallax tutorial 2018

Coryn Bailer-Jones, MPIA Heidelberg (https://mpia.de/homes/calj)

This tutorial concerns the use of inferring distances and velocities from parallaxes and proper motions using Bayesian inference. It uses R codes in jupyter notebooks, together with simulated data or data from Gaia-DR1 (TGAS). The tutorial is divided into three parts, each with its own notebook:

1. Inference of distance to a single source using its parallax. Includes: various priors, simple hierarchical model. 
2. Inference of distance to and size of a cluster using parallax and positions of its members. Includes: accommodating to correlated measurements. See resources/cluster_inference.pdf for details.
3. Inference of distance to and 2D tangential velocity on the sky of a single source using its parallax and proper motion. Includes: use of MCMC. See resources/3D_astrometry_inference.pdf for details. 

Resources amnd further information:
* [Cluster distance inference] (resources/cluster_inference.pdf)
* [Joint inference from parallax and proper motions (CBJ-081)] (resources/3D_astrometry_inference.pdf)
* [Gaia Data Release 1](http://adsabs.harvard.edu/abs/2017A%26A...601A..19G)
* [parallax paper 1](http://adsabs.harvard.edu/abs/2015PASP..127..994B)
* [parallax paper 2](http://adsabs.harvard.edu/abs/2016ApJ...832..137A)
* [parallax paper 3](http://adsabs.harvard.edu/abs/2016ApJ...833..119A)
* [Practical Bayesian Inference](http://www.cambridge.org/de/academic/subjects/physics/mathematical-methods/practical-bayesian-inference-primer-physical-scientists?format=PB)

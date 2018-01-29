# Parallax tutorial 2018

Coryn Bailer-Jones, MPIA Heidelberg (https://mpia.de/homes/calj)

This tutorial concerns the use of inferring distances and velocities from parallaxes and proper motions using Bayesian inference. It uses R codes in jupyter notebooks, together with simulated data or data from Gaia-DR1 (TGAS). This tutorial is divided into three parts, each with its own notebook:

1. Inference of distance to a single source using its parallax
2. Inference of distance to and size of a cluster using parallax and positions of its members. This includes accommodating to correlated measurements
3. Inference of distance to and 2D tangential velocity on the sky of a single source using its parallax and proper motion 

The PDF file cluster_inference.pdf explains more precisely the approach for inferring cluster distance and size/shape.

There are three python notebooks with exercises which illustrate the statistical methods covered in the lectures.
In your root directory put the notebooks (*.ipynb) into a directory called "notebook/", the R codes (*.R) into a directory called "Rcode/", and the data (*.csv) in a directory called "data/". The root directory is then defined in the first line of each notebook.

Resources:
* [Cluster distance inference] (cluster_inference.pdf)
* [Gaia Data Release 1](http://adsabs.harvard.edu/abs/2017A%26A...601A..19G)
* [parallax paper 1](http://adsabs.harvard.edu/abs/2015PASP..127..994B)
* [parallax paper 2](http://adsabs.harvard.edu/abs/2016ApJ...832..137A)
* [parallax paper 3](http://adsabs.harvard.edu/abs/2016ApJ...833..119A)
* [Practical Bayesian Inference](http://www.cambridge.org/de/academic/subjects/physics/mathematical-methods/practical-bayesian-inference-primer-physical-scientists?format=PB)

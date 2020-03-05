# thsiao
[![Build Status](https://travis-ci.com/tXiao95/thsiao.svg?branch=master)](https://travis-ci.com/tXiao95/thsiao)

Personal convenience functions for command line and numerical tools. Also for package demonstrative purposes.

## Modeling tools
While there are [packages for simulating responses from a distribution corresponding to a fitted model object](https://cran.r-project.org/web/packages/arm/arm.pdf), and functions that
predict responses given new data, there are no methods that simulate responses for new data. Such simulation is important in 
Monte Carlo methods. In the completely linear case, standard matrix multiplication is simple but code becomes more complicated when splines like `splines::bs` or 
other transformations are added into the formula. `simNew` is a `S3` method currently defined for `lm` and `glm` objects that simulates responses
on new data not originally used to fit the model object. For example, for a `lm` object, `simNew` will sample from the multivariate normal distribution 
defined by the coefficients of the model, and apply those samples to the new data, returning a `m x n` matrix, where `m` is the number of rows 
in the new dataset, and `n` is the number of simulations requested.

## MATLAB gallery of test matrices
Bringing the MATLAB gallery of test matrices into R. https://www.mathworks.com/help/matlab/ref/gallery.html

## Command line and SGE cluster tools 
SGE cluster tools like easier profiling, monitoring, and manipulating of jobs from R. Also command line interfacing from R.

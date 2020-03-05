# thsiao
[![Build Status](https://travis-ci.com/tXiao95/thsiao.svg?branch=master)](https://travis-ci.com/tXiao95/thsiao)

Personal convenience functions for command line and numerical tools. Also for package demonstrative purposes.

## Modeling tools
While there are packages for simulating responses from a distribution corresponding to a fitted model object, and functions that
predict responses given new data, there are no methods that simulate responses for new data. For example, for an \code{lm} object, 
\code{simulateDraws} will sample from the multivariate normal distribution defined by the coefficients of the model, and apply those
samples to the new data generating a $m \times n$ matrix, where $m$ is the number of rows in the new dataset, and $n$ is the number of 
simulations requested.

## MATLAB gallery of test matrices
Bringing the MATLAB gallery of test matrices into R. https://www.mathworks.com/help/matlab/ref/gallery.html

## Command line and SGE cluster tools 

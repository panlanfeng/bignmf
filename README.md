bignmf
======

This is a R package intending to solve Nonnegative Matrix Factorization efficiently. 
The algorithm bases on alternating least squares while solving nonnegative constrained regression via coordinate descent. 
The coordinate descent method is modified from glmnet (See Regularization Paths for Generalized Linear Models via Coordinate Descent, by Friedman et al. 2010). 

The inner part of the algorithm is implemented in C++. 
This package currently tested only on Linux. 


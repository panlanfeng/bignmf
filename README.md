bignmf
======

This R package intends to solve Nonnegative Matrix Factorization efficiently. 
The algorithm bases on alternating least squares while solving nonnegative constrained regression via coordinate descent. 
The coordinate descent method is modified from glmnet (See __Regularization Paths for Generalized Linear Models via Coordinate Descent__, by J Friedman et al. 2010). 

The inner part of the algorithm is implemented in C++. 
This package currently tested only on Linux. 

##How to install

In Linux teminal, go to some directory:  

      	  cd  ~/mypackages/ 

Download it with:  

      	 git clone https://github.com/panlanfeng/bignmf.git

Build with:  

    	  R CMD build bignmf

Then in R:  

    	  install.packages("~/mypackages/bignmf_0.1.tar.gz", repos = NULL, type = "source")   
	  library(bignmf)  
	  
	  
	


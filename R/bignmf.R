##V: the matrix to be factorized.
##r: the rank of resulting matrices.
##initial: can be "H.kmeans", "W.kmeans", or the initial W or H matrix. "H.kmeans" means initialize H with kmeans(V,r)$centers, while "W.kmeans" means initialized W with kmeans centers. The default is "H.kmeans".
##max.iteration: the number of iterations allowed.
## stop.condition: the function compares the norm of projected gradient matrix in the k-th iteration and the norm of gradient matrix after the first iteration. If the former one is less than the latter multiplying stop.condition, iteration stops .

##Detail: The nonnegative matrix factorization tries to find nonnegative matrices W and H, so that V \approx WH. Using sum of squares loss function, the problem is to solve \min_{W\ge0, H\ge0} f(V - WH). bignmf finds W minimizing f given H and then finds H give W, i.e. alternating least squares. The function treats nonnegative constrained regression as a special L1 regression and solves it via coordinate descent method. 


##value: the function returns a list of length 3.
##W : the resulting nonnegative matrix W.
##H : the resulting nonnegative matrix H.
##iterations : number of iterations.

##Example:
# v_mat <- matrix(rexp(60000,2), 200, 300)
# system.time(re <- bignmf(v_mat, 5))
# re$iterations
# 
# v_mat <- matrix(rexp(6000000,2), 2000, 3000)
# v_mat[v_mat < quantile(v_mat, .1)] <- 0
# system.time(re <- bignmf(v_mat, 20))
# re$iterations


bignmf <- function(V, r=5, max.iteration=200, stop.condition=1e-4){
  
  V <- as.matrix(V)
  if(storage.mode(V)!="double"){
    storage.mode(V) <- "double"
  }
  
  nm <- dim(V)
  n <- nm[1]
  m <- nm[2]
  
  W <- abs(matrix(rnorm(r * n), n, r))
  H <- abs(matrix(rnorm(r * m), r, m))
  
  wh <- .Call("whupdate",V, W, H, as.integer(max.iteration), stop.condition)
  
  if(wh$iterations == max.iteration)
    warning("Iteration doesn't converge!")
  
  return(wh)
}


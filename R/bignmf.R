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


bignmf <- function(V, r, initial="H.random", max.iteration=200, stop.condition=1e-3){
  
  V <- as.matrix(V)
  if(storage.mode(V)!="double"){
    storage.mode(V) <- "double"
  }
  
  nm <- dim(V)
  n <- nm[1]
  m <- nm[2]
  eps <- .Machine$double.eps
  
  if(is.matrix(initial)){
    initial.dim <- dim(initial)
    if(initial.dim[1] == n & initial.dim[2] == r){
      W <- initial
      H <- .Call("hupdate", V, W)
    }    else if(initial.dim[1] == r & initial.dim[2] == m){
      H <- initial
      #W <- .Call("wupdate", V, H)
    }    else{
      stop("The initial value is of wrong dimension!")
    }
  } else if (class(initial) == "character") {
    if(initial == "H.kmeans"){
      H <- as.matrix(kmeans(V, r)$centers)
      #        W <- matrix(0, n, r)
    } else if(initial == "W.kmeans"){
      W<- as.matrix(t(kmeans(t(V), r)$centers))
      H <- .Call("hupdate", V, W)
    } else if(initial == "H.random"){
      H <- matrix(rexp(r * m, 1), r, m)
    }  else if(initial == "W.random"){
      W <- matrix(rexp(r * n, 1), n, r)
      H <- .Call("hupdate", V, W)
    } else{
      stop("Please give the right initial value!")
    }
    
  }else{
    stop("initial should be a matrix or a character!")
  }
  
  
  
  
  for(i in 1:max.iteration){
    
    #Update W  
    W <- .Call("wupdate",V, H)
    # Update H
    H <- .Call("hupdate",V, W)
    
    #Keep W and H in a similar scale
    dk <- sqrt(sqrt(rowSums(H ^ 2)) / (sqrt(colSums(W ^ 2)) + .Machine$double.eps))
    W <- W %*% diag(dk)
    H <- diag(1 / (dk + .Machine$double.eps)) %*% H
    
    
    #check whether the iteration  converge or not. 
    #This is the criterion is based on Projected Gradient(Lin,2007).
    if(i == 1){
      res <- W %*% H - V
      grad.w1h1 <- sum(tcrossprod(res, H) ^ 2) + sum(crossprod(W, res) ^ 2)
    }else if(i %% 5 == 0) {
      res <- W %*% H - V
      grad.w <- tcrossprod(res, H)
      # if W_ij=0, grad.w_ij=min(0, grad.w_ij)
      grad.w[W == 0] <- grad.w[W == 0] * (grad.w[W == 0] < 0)
      grad.h <- crossprod(W, res)
      grad.h[H == 0] <- grad.h[H == 0] * (grad.h[H == 0] < 0)
      grad.wh <- sum(grad.w ^ 2) + sum(grad.h ^ 2)
      
      if(grad.wh < grad.w1h1 * stop.condition){
        break
      } 
    }
    
    
    if(i == max.iteration)
      warning("Iteration doesn't converge!")
  }
  
  return(list(W=W, H=H, iterations=i))
}


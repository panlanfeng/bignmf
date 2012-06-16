\name{bignmfsp}
\alias{bignmfsp}
\title{
Solving the sparse nonnegative matrix factorization via alternating least square.
}
\description{
Solving the sparse nonnegative matrix factorization via alternating least square.
}
\usage{
bignmfsp(V, r, initial = "H.random", max.iteration = 200, stop.condition = 0.001)
}


\arguments{
  \item{V}{
the matrix to be factorized. Should be of class "dgCMatrix" from Matrix Package.
}
  \item{r}{
the rank of resulting matrices.
}
  \item{initial}{
can be \code{"H.kmeans"}, \code{"W.kmeans"}, \code{"H.random"}, \code{"W.random"}, or a matrix giving initial W or H. 
\code{"H.kmeans"} means initialize H with \code{kmeans(V,r)$centers}, while 
\code{"W.kmeans"} means 
initialized W with kmeans centers. \code{"H.random"} means initialize H with random numbers from an exponential distribution, while 
\code{"W.random"} means 
initialized W randomly. The default is \code{"H.random"}.
}
  \item{max.iteration}{
the number of iterations allowed.
}
  \item{stop.condition}{
the function compares the norm of projected gradient matrix in the k-th iteration 
and the norm of gradient matrix after the first iteration. If the former one is less
than the latter multiplying stop.condition, iteration stops .
}
}
\details{
The nonnegative matrix factorization tries to find nonnegative matrices W and H, so that 
V \eqn{\approx} WH. Using sum of squares loss function, the problem is to solve 
\eqn{\min_{W\ge0, H\ge0} f(V - WH)}.
bignmf finds W minimizing f given H and then finds H give W, i.e. alternating least squares. 
The function treats nonnegative constrained regression as a special L1 regression and solves it 
via coordinate descent method. 

}
\value{
The function returns a list of length 3.
\item{W}{ the resulting nonnegative matrix W.}
\item{H}{ the resulting nonnegative matrix H.}
\item{iterations}{ number of iterations.}
}


% \seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{
v_mat <- matrix(rexp(60000,2), 200, 300)
system.time(re <- bignmf(v_mat, 5))
re$iterations
 
v_mat <- matrix(rexp(6000000,2), 2000, 3000)
v_mat[v_mat < quantile(v_mat, .1)] <- 0
system.time(re <- bignmf(v_mat, 20))
re$iterations
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{optimize}


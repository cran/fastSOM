\name{soi_avg_exact}
\Rdversion{1.1}
\alias{soi_avg_exact}
\title{Exact Calculation of Average, Minimal, and Maximal Spillover Index}
\description{
  Calculates the Average, Minimal, and Maximal Spillover Index exactly.
}
\usage{soi_avg_exact(Sigma, A, ncores = 1)}
\arguments{
  \item{Sigma}{ 
    Either a covariance matrix or a list thereof.
  }
  \item{A}{ 
    Either a 3-dimensional array with A[,,h] being MA coefficient matrices of the same dimension as \code{Sigma} or a list thereof.
  }
  \item{ncores}{ 
	Number of cores. Missing ncores or \code{ncores=1} means no parallelization (just one core is used). 
    \code{ncores=0} means automatic detection of the number of available cores. 
    Any other integer determines the maximal number of cores to be used.}
}
\details{
  The spillover index introduced by Diebold and Yilmaz (2009) (see References) depends on the ordering of the model variables.
  While \code{\link{soi_avg_est}} 
  provides an algorithm to estimate average, minimum, and maximum of the spillover index over all permutations, 
  \code{soi_avg_est} calculates these quantities exactly. Notice, however, that for large dimensions \eqn{N}, this might be quite
  time- as well as memory-consuming.
  If only the exact average of the spillover index is wanted, \code{soi_from_sot(sot_avg_exact(Sigma,A,ncores)$Average)} should be used.
  
  The typical application of the 'list' version of \code{soi_avg_exact} is a rolling windows approach when \code{Sigma} and \code{A} are lists representing the corresponding quantities at different points in time
  (rolling windows). 
}
\value{
  The 'single' version returns a list containing the exact average, minimal, and maximal spillover index as well as permutations that generated the minimal and maximal value. 
  The 'list' version returns a list consisting of three vectors (the average, minimal, and maximal spillover index values) and two matrices (the columns of which are the permutations generating the minima and maxima). 
}
\references{
  [1] Diebold, F. X. and Yilmaz, K. (2009):  \href{http://onlinelibrary.wiley.com/doi/10.1111/j.1468-0297.2008.02208.x/pdf}{Measuring financial asset return and volatitliy spillovers, 
                                             with application to global equity markets}, 
                                             Economic Journal 199(534): 158-171.

  [2] Kloessner, S. and Wagner, S. (2012):  \href{http://onlinelibrary.wiley.com/doi/10.1002/jae.2366/pdf}{Exploring All VAR Orderings for Calculating Spillovers? Yes, We Can! - 
											A Note on Diebold and Yilmaz (2009)}, Journal of Applied Econometrics 29(1): 172-179
}
\author{
  Stefan Kloessner (\email{S.Kloessner@mx.uni-saarland.de}), \cr
  with contributions by Sven Wagner (\email{sven.wagner@mx.uni-saarland.de})
}
\seealso{ \code{\link{fastSOM-package}}, \code{\link{soi_avg_est}} }
\examples{
# generate randomly positive definite matrix Sigma of dimension N 
N <- 10
Sigma <- crossprod(matrix(rnorm(N*N),nrow=N)) 
# generate randomly coefficient matrices
H <- 10 
A <- array(rnorm(N*N*H),dim=c(N,N,H)) 
# calculate the exact average, minimal, 
# and maximal spillover index and determine the corresponding ordering 
# of the model variables
soi_avg_exact(Sigma, A) 
}
\keyword{exact average of spillover index}
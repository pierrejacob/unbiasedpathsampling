# Couplings of Normal distributions, adapted from the 'unbiasedpoisson' package.
# All functions return a list with 'xy' (the two draws) and 'identical' (TRUE if x = y).

#'@rdname rmvnorm_max
#'@title Maximal coupling of two multivariate Normal distributions
#'@description Sample from maximal coupling of two multivariate Normal distributions,
#'specified through their means and covariance matrices. See \code{\link{rmvnorm_max_chol}}
#'for a version using Cholesky factors and their inverses.
#'See \code{\link{rmvnorm_reflectionmax}}
#'for a reflection-maximal coupling, in the case Sigma1=Sigma2.
#'@param mu1 First mean
#'@param mu2 Second mean
#'@param Sigma1 First covariance matrix
#'@param Sigma2 Second covariance matrix
#'@return A list containing 'xy', a matrix with 2 columns (one for each draw),
#' and a boolean indicator 'identical' indicating whether the two draws
#' are identical.
#'@examples
#' p <- 3
#' mu1 <- rep(0, p)
#' mu2 <- rep(1, p)
#' Sigma1 <- diag(0.4, p, p)
#' Sigma1[1,2] <- Sigma1[2,1] <- 0.2
#' Sigma2 <- diag(1.4, p, p)
#' Sigma2[1,2] <- Sigma2[2,1] <- -0.5
#' rmvnorm_max(mu1, mu2, Sigma1, Sigma2)
#'@export
rmvnorm_max <- function(mu1, mu2, Sigma1, Sigma2){
  return(rmvnorm_max_coupling_(mu1, mu2, Sigma1, Sigma2))
}

#'@rdname rmvnorm_max_chol
#'@title Maximal coupling of two multivariate Normal distributions
#'@description Sample from maximal coupling of two multivariate Normal distributions,
#'specified through their means, the Cholesky factors of their covariance matrices,
#'and the inverse of the Cholesky factors of the covariance matrices.
#'@param mu1 First mean
#'@param mu2 Second mean
#'@param Cholesky1 First Cholesky factor, e.g. obtained with \code{\link[base]{chol}}
#'@param Cholesky2 Second Cholesky factor
#'@param Cholesky_inverse1 First inverse of Cholesky factor, e.g. obtained with \code{solve(chol(Sigma))}
#'@param Cholesky_inverse2 Second inverse of Cholesky factor
#'@return A list containing 'xy', a matrix with 2 columns (one for each draw),
#' and a boolean indicator 'identical' indicating whether the two draws
#' are identical.
#'@examples
#' p <- 3
#' mu1 <- rep(0, p)
#' mu2 <- rep(1, p)
#' Sigma1 <- diag(0.4, p, p)
#' Sigma1[1,2] <- Sigma1[2,1] <- 0.2
#' Sigma2 <- diag(1.4, p, p)
#' Sigma2[1,2] <- Sigma2[2,1] <- -0.5
#' Sigma1_chol <- chol(Sigma1)
#' Sigma2_chol <- chol(Sigma2)
#' Sigma1_chol_inv <- solve(Sigma1_chol)
#' Sigma2_chol_inv <- solve(Sigma2_chol)
#' rmvnorm_max_chol(mu1, mu2, Sigma1_chol, Sigma2_chol, Sigma1_chol_inv, Sigma2_chol_inv)
#'@export
rmvnorm_max_chol <- function(mu1, mu2, Cholesky1, Cholesky2, Cholesky_inverse1, Cholesky_inverse2){
  # we need Cholesky <- chol(Sigma), not the transpose
  valid_Cholesky1 <- all(Cholesky1[lower.tri(Cholesky1)]==0)
  valid_Cholesky2 <- all(Cholesky2[lower.tri(Cholesky2)]==0)
  stopifnot(valid_Cholesky1, valid_Cholesky2)
  return(rmvnorm_max_coupling_cholesky_(mu1, mu2, Cholesky1, Cholesky2, Cholesky_inverse1, Cholesky_inverse2))
}

#'@rdname rmvnorm_reflectionmax
#'@title Reflection-maximal coupling of two multivariate Normal distributions
#'@description Sample from reflection-maximal coupling of two multivariate Normal distributions,
#'specified through their means, with the same covariance matrix, specified
#'through its Cholesky factor and inverse of Cholesky factor.
#'
#'The idea is that a multivariate Normal is drawn around the first mean (mu1),
#'and then reflected with respect to a hyperplane orthogonal to the direction between mu1 and mu2.
#'
#'For univariate Normal distributions, see \code{\link{rnorm_reflectionmax}};
#'for diagonal covariance matrices, see \code{\link{rmvnorm_reflectionmax_diag}}.
#'
#'@param mu1 First mean
#'@param mu2 Second mean
#'@param Cholesky Cholesky factor, e.g. obtained with \code{\link[base]{chol}}
#'@param Cholesky_inverse Inverse of Cholesky factor, e.g. obtained with \code{solve(chol(Sigma))}
#'@return A list containing 'xy', a matrix with 2 columns (one for each draw),
#' and a boolean indicator 'identical' indicating whether the two draws
#' are identical.
#'@examples
#' p <- 3
#' mu1 <- rep(0, p)
#' mu2 <- rep(1, p)
#' Sigma <- diag(0.4, p, p)
#' Sigma[1,2] <- Sigma[2,1] <- 0.2
#' Sigma_chol <- chol(Sigma)
#' Sigma_chol_inv <- solve(Sigma_chol)
#' rmvnorm_reflectionmax(mu1, mu2, Sigma_chol, Sigma_chol_inv)
#'@export
rmvnorm_reflectionmax <- function(mu1, mu2, Cholesky, Cholesky_inverse){
  d_ <- dim(Cholesky)
  if (is.null(d_) || (d_[1] == 1)){
    stop("function 'rmvnorm_reflectionmax' is meant for multivariate Normals; for univariate Normals, use 'rnorm_reflectionmax'")
  }
  return(rmvnorm_reflection_max_coupling_(mu1, mu2, Cholesky, Cholesky_inverse))
}

#'@rdname rmvnorm_reflectionmax_diag
#'@title Reflection-maximal coupling of two multivariate Normal distributions with diagonal covariance
#'@description Sample from reflection-maximal coupling of two multivariate Normal distributions,
#'specified through their means, with the same diagonal covariance matrix, specified
#'through the vector of standard deviations (the diagonal of its Cholesky factor).
#'
#'@param mu1 First mean
#'@param mu2 Second mean
#'@param Cholesky vector of standard deviations, i.e. square roots of the diagonal of the covariance matrix
#'@return A list containing 'xy', a matrix with 2 columns (one for each draw),
#' and a boolean indicator 'identical' indicating whether the two draws
#' are identical.
#'@examples
#' p <- 3
#' mu1 <- rep(0, p)
#' mu2 <- rep(1, p)
#' Sigma <- rep(0.4^2, p)
#' Sigma_chol <- sqrt(Sigma)
#' rmvnorm_reflectionmax_diag(mu1, mu2, Sigma_chol)
#'@export
rmvnorm_reflectionmax_diag <- function(mu1, mu2, Cholesky){
  return(rmvnorm_reflection_max_coupling_diag_(mu1, mu2, Cholesky))
}

#'@rdname rnorm_reflectionmax
#'@title Reflection-maximal coupling of two univariate Normal distributions
#'@description Sample from reflection-maximal coupling of two univariate Normal distributions,
#'specified through their means, with common standard deviation.
#'See \code{\link{rmvnorm_reflectionmax}} for a multivariate version.
#'@param mu1 First mean
#'@param mu2 Second mean
#'@param sigma Common standard deviation
#'@return Returns a list with
#' \itemize{
#'
#' \item "xy": the pair of samples \eqn{(x,y)}
#'
#' \item "identical": TRUE if \eqn{x = y}, FALSE otherwise
#' }
#'@examples
#' rnorm_reflectionmax(0, 1, 0.5)
#'@export
rnorm_reflectionmax <- function(mu1, mu2, sigma){
  reflmax_xy <- c(0,0)
  # draw std normal first
  xdot <- rnorm(1)
  # this follows the notation of Bou Rabee et al, 2018, roughly
  z <- (mu1 - mu2) / sigma
  normz <- abs(z)
  ident_ <- FALSE
  if (normz < 1e-15){
    ident_ <- TRUE
  } else {
    e <- z / normz
    utilde <- runif(1, 0, 1)
    ident_ <- (log(utilde) < (dnorm(xdot + z, 0, 1, log = TRUE) - dnorm(xdot, log = TRUE)))
  }
  reflmax_xy[1] <- mu1 + sigma * xdot
  if (ident_){
    # x = y exactly (avoids floating-point discrepancies)
    reflmax_xy[2] <- reflmax_xy[1]
  } else {
    ydot <- xdot - 2 * (e * xdot) * e
    reflmax_xy[2] <- mu2 + sigma * ydot
  }
  return(list(xy = reflmax_xy, identical = ident_))
}

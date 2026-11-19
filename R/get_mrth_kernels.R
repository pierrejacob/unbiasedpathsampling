#'@rdname get_mrth_kernels
#'@title Get random walk Metropolis-Rosenbluth-Teller-Hastings kernels
#'@description This function takes a target (specified through its log-pdf)
#' and a covariance matrix for a Normal random walk proposal, and returns a list containing the keys
#' \code{single_kernel}, \code{coupled_kernel} corresponding to marginal
#' and coupled MRTH kernels. The target can be univariate or multivariate.
#'
#' The coupling is done by reflection-maximal coupling of the proposals,
#' and a common uniform variable for the accept/reject step. For reflection-maximal
#' couplings, see \code{\link{rnorm_reflectionmax}}, \code{\link{rmvnorm_reflectionmax}}
#' and \code{\link{rmvnorm_reflectionmax_diag}}.
#'
#' The returned kernels can then be used in the functions \code{\link{coupled_chains}}
#' or \code{\link{unbiased_estimator}}.
#'
#'@param logtarget function taking a vector as input and returning target log-density evaluation
#'@param Sigma_proposal covariance of the Normal random walk proposal. Either
#' a d x d covariance matrix, or a vector of length d containing the diagonal of a diagonal
#' covariance matrix, or a single number (the variance) for univariate targets.
#'@return A list containing the keys
#' \code{single_kernel}, \code{coupled_kernel}.
#'
#' \code{single_kernel(chain_state, current_pdf)} returns a list with keys
#' \code{chain_state}, \code{current_pdf} and \code{accept}.
#'
#' \code{coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2)} returns a list with keys
#' \code{chain_state1}, \code{chain_state2}, \code{current_pdf1}, \code{current_pdf2}
#' and \code{identical} (TRUE if the two chains are equal after the transition).
#'@examples
#' Sigma_target <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' target <- get_mvnormal(2, c(0, 0), Sigma_target)
#' kernels <- get_mrth_kernels(target$logtarget, Sigma_target)
#' rinit <- function(){
#'   chain_state <- rnorm(2)
#'   list(chain_state = chain_state, current_pdf = target$logtarget(chain_state))
#' }
#' unbiased_estimator(kernels$single_kernel, kernels$coupled_kernel, rinit, k = 10, m = 50)
#'@export
get_mrth_kernels <- function(logtarget, Sigma_proposal){
  ## work out the type of proposal covariance
  if (is.matrix(Sigma_proposal)){
    dimension <- nrow(Sigma_proposal)
    stopifnot(ncol(Sigma_proposal) == dimension)
    if (dimension == 1){
      type <- "univariate"
      proposal_sd <- sqrt(Sigma_proposal[1,1])
    } else if (all(Sigma_proposal[row(Sigma_proposal) != col(Sigma_proposal)] == 0)){
      # diagonal matrix: use the faster diagonal version
      type <- "diagonal"
      proposal_sd <- sqrt(diag(Sigma_proposal))
    } else {
      type <- "full"
      Sigma_proposal_chol <- chol(Sigma_proposal)
      Sigma_proposal_chol_inv <- solve(Sigma_proposal_chol)
    }
  } else {
    dimension <- length(Sigma_proposal)
    if (dimension == 1){
      type <- "univariate"
    } else {
      type <- "diagonal"
    }
    proposal_sd <- sqrt(Sigma_proposal)
  }
  stopifnot(all(is.finite(Sigma_proposal)))

  ## sample from the proposal given the current state
  propose <- switch(type,
                    univariate = function(chain_state) rnorm(1, chain_state, proposal_sd),
                    diagonal = function(chain_state) chain_state + proposal_sd * rnorm(dimension),
                    full = function(chain_state) fast_rmvnorm_chol(1, chain_state, Sigma_proposal_chol)[1,])
  ## sample from the reflection-maximal coupling of the proposals given two states
  coupled_propose <- switch(type,
                            univariate = function(chain_state1, chain_state2){
                              rnorm_reflectionmax(chain_state1, chain_state2, proposal_sd)
                            },
                            diagonal = function(chain_state1, chain_state2){
                              rmvnorm_reflectionmax_diag(chain_state1, chain_state2, proposal_sd)
                            },
                            full = function(chain_state1, chain_state2){
                              rmvnorm_reflectionmax(chain_state1, chain_state2, Sigma_proposal_chol, Sigma_proposal_chol_inv)
                            })

  # single kernel
  single_kernel <- function(chain_state, current_pdf){
    proposal_value <- propose(chain_state)
    proposal_pdf <- logtarget(proposal_value)
    accept <- is.finite(proposal_pdf) && (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(list(chain_state = proposal_value, current_pdf = proposal_pdf, accept = accept))
    } else {
      return(list(chain_state = chain_state, current_pdf = current_pdf, accept = accept))
    }
  }

  # coupled kernel
  coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2){
    proposal_value <- coupled_propose(chain_state1, chain_state2)
    if (type == "univariate"){
      proposal1 <- proposal_value$xy[1]
      proposal2 <- proposal_value$xy[2]
    } else {
      proposal1 <- proposal_value$xy[,1]
      proposal2 <- proposal_value$xy[,2]
    }
    proposal_pdf1 <- logtarget(proposal1)
    if (proposal_value$identical){
      # no need to evaluate the target twice
      proposal2 <- proposal1
      proposal_pdf2 <- proposal_pdf1
    } else {
      proposal_pdf2 <- logtarget(proposal2)
    }
    # common uniform for the accept/reject steps
    logu <- log(runif(1))
    accept1 <- FALSE
    accept2 <- FALSE
    if (is.finite(proposal_pdf1)){
      accept1 <- (logu < (proposal_pdf1 - current_pdf1))
    }
    if (is.finite(proposal_pdf2)){
      accept2 <- (logu < (proposal_pdf2 - current_pdf2))
    }
    if (accept1){
      chain_state1 <- proposal1
      current_pdf1 <- proposal_pdf1
    }
    if (accept2){
      chain_state2 <- proposal2
      current_pdf2 <- proposal_pdf2
    }
    identical_ <- proposal_value$identical && accept1 && accept2
    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2,
                current_pdf1 = current_pdf1, current_pdf2 = current_pdf2,
                identical = identical_))
  }
  return(list(single_kernel = single_kernel, coupled_kernel = coupled_kernel))
}

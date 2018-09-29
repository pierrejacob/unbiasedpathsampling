# load package
library(unbiasedpathsampling)
# load packages for parallel computation
library(doParallel)
library(doRNG)
# load package for data manipulation
library(dplyr)
# load packages for plotting
library(ggplot2)
library(latex2exp)
# load graphical theme
setmytheme()
# register processing cores
registerDoParallel(cores = detectCores()-2)
# clean environment
rm(list = ls())
# set seed
set.seed(21)

p <- 7 # number of covariates
theta_star <- (0:6)/10 # data-generating parameter
expit <- function(z) 1 / (1 + exp(-z)) # expit function
n <- 1000 # number of observations
# generate X and Y from the model
X <- cbind(rep(1, n), scale(matrix(rnorm(n*(p-1)), ncol = p-1)))
logitprobs <- apply(X = X, MARGIN = 1, FUN = function(row) sum(row * theta_star))
Y <- rbinom(n, 1, expit(logitprobs))
b <- matrix(0, nrow = p, ncol = 1)
B <- diag(10, p, p)

## Check with Importance sampling
library(Rcpp)
cppFunction('
            double loglikelihood_cpp(const NumericVector & Y, const NumericMatrix & X, const NumericVector & theta){
            double logpr = 0.;
            int p = theta.size();
            int ny = X.rows();
            for (int i = 0; i < ny; i++){
            double thetaX = 0;
            for (int j = 0; j < p; j++){
            thetaX = thetaX + theta(j) * X(i,j);
            }
            logpr = logpr - log(1+exp(-thetaX)) + (1-Y(i)) * (-thetaX);
            }
            return(logpr);
            }')

logisticmle <- optim(par = theta_star, fn = function(theta) - loglikelihood_cpp(Y, X, theta))
theta_mle <- logisticmle$par
library(numDeriv)
hess <- hessian(func = function(theta) - loglikelihood_cpp(Y, X, theta), x = theta_mle)
asymptotic_var <- solve(hess)

## estimate cst with IS
NIS <- 1e4
samplesIS <- fast_rmvnorm(NIS, theta_mle, asymptotic_var)
postdensities <- apply(samplesIS, 1, function(x) loglikelihood_cpp(Y, X, x) + fast_dmvnorm(t(x), b, B))
logweights <- postdensities - fast_dmvnorm(samplesIS, theta_mle, asymptotic_var)
maxweights <- max(logweights)
weights <- exp(logweights-maxweights)
cst <- log(mean(weights)) + maxweights
normweights <- weights / sum(weights)
1/(sum(normweights^2)) # ESS
# log Z_1
cst
# log Z_1 - log Z_0 where Z_0 = 2^{-n}
cst + n * log(2)

##
## Now path from Normal approx to posterior

target0 <- function(x) fast_dmvnorm(t(x), theta_mle, asymptotic_var)
target1 <- function(x) loglikelihood_cpp(Y, X, x) + fast_dmvnorm(t(x), b, B)
target_theta <- function(x, theta) (1 - theta) * target0(x) + theta * target1(x)
loggammaratio <- function(x) target1(x) - target0(x)

dimension <- 7

rinit <- function(theta){
  chain_state <- fast_rmvnorm(1, mean = theta_mle, covariance = asymptotic_var)[1,]
  current_pdf <- target_theta(chain_state, theta)
  return(list(chain_state = chain_state, current_pdf = current_pdf))
}

# Metropolis hastings for a certain theta between 0 and 1
Sigma_q <- asymptotic_var/dimension
Sigma_chol <- chol(Sigma_q)
Sigma_chol_inv <- solve(Sigma_chol)

single_kernel <- function(chain_state, current_pdf, theta){
  proposal_value <- chain_state + fast_rmvnorm_chol(1, rep(0, dimension), Sigma_chol)[1,]
  proposal_pdf <- target_theta(proposal_value, theta)
  accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
  if (accept){
    return(list(chain_state = proposal_value, current_pdf = proposal_pdf))
  } else {
    return(list(chain_state = chain_state, current_pdf = current_pdf))
  }
}

# Markov kernel of the coupled chain
coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, theta){
  cproposals_ <- unbiasedpathsampling:::rnorm_reflectionmax_(chain_state1, chain_state2, Sigma_chol, Sigma_chol_inv)
  proposal1 <- cproposals_$xy[,1]
  proposal2 <- cproposals_$xy[,2]
  proposal_pdf1 <- target_theta(proposal1, theta)
  proposal_pdf2 <- target_theta(proposal2, theta)
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
  return(list(chain_state1 = chain_state1, chain_state2 = chain_state2,
              current_pdf1 = current_pdf1, current_pdf2 = current_pdf2))
}

theta <- 0
ri <- function() rinit(theta)
sk <- function(x, f) single_kernel(x, f, theta)
ck <- function(x1, x2, f1, f2) coupled_kernel(x1, x2, f1, f2, theta)
nrep <- 1000
meetings <- foreach(i = 1:nrep, .combine = c) %dorng% {
  unbiased_estimator(sk, ck, ri)$meetingtime
}
k <- as.numeric(floor(quantile(meetings, probs = 0.99)))
m <- 5*k

uestimator_giventheta <- function(theta){
  ri <- function() rinit(theta)
  sk <- function(x, f) single_kernel(x, f, theta)
  ck <- function(x1, x2, f1, f2) coupled_kernel(x1, x2, f1, f2, theta)
  ue_ <- unbiased_estimator(sk, ck, ri, h = loggammaratio, k = k, m = m)
  return(ue_)
}

UPS <- function(rp, dp, n){
  thetas <- rp(n)
  results.df <- foreach(i = 1:n, .combine = rbind) %dorng% {
    theta <- thetas[i]
    unbiased_result <- uestimator_giventheta(theta)
    uestimator <- unbiased_result$uestimator[[1]] / dp(theta) # unbiased estimator
    meetingtime <- unbiased_result$meetingtime # meeting time tau
    iteration <- unbiased_result$iteration # total number of iteration (max(m,tau))
    cost <- 2*meetingtime + max(1, iteration + 1 - meetingtime)
    data.frame(i = 1, theta = theta, uestimator = uestimator, meetingtime = meetingtime, iteration = iteration, cost = cost)
  }
  return(results.df)
}

nrep <- 100
results.optimal.df <- UPS(runif, dunif, nrep)
g <- ggplot(results.optimal.df, aes(x = theta, y = uestimator)) + geom_point(alpha = 0.25)
g <- g + xlab(expression(theta)) + ylab(TeX("$\\hat{E}(\\theta)/q(\\theta)$"))
g
# ggsave(filename = "ups.logistic.otherpath.png", plot = g, width = 5, height = 5)

## Confidence intervals
logZ1overZ0_hat_opt <- mean(results.optimal.df$uestimator)
logZ1overZ0_variance_opt <- var(results.optimal.df$uestimator)
logZ1overZ0_cost_opt <- mean(results.optimal.df$cost)
logZ1overZ0_error_opt <- sqrt(logZ1overZ0_variance_opt/nrow(results.optimal.df))
cat("[", logZ1overZ0_hat_opt-1.96*logZ1overZ0_error_opt, ",", logZ1overZ0_hat_opt+1.96*logZ1overZ0_error_opt, "]\n")
cst

cat("[", n*log(2) + logZ1overZ0_hat_opt-1.96*logZ1overZ0_error_opt, ",", n*log(2) + logZ1overZ0_hat_opt+1.96*logZ1overZ0_error_opt, "]\n")
cst + n*log(2)

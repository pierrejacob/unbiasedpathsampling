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
library(MASS)
data("mammals")
gdata <- qplot(mammals$body, mammals$brain) + scale_x_log10() + scale_y_log10()
gdata <- gdata + xlab("body") + ylab("brain")
gdata
# ggsave(filename = "ups.cv.mammals.data.pdf", plot = gdata, width = 6, height = 5)

Xfull <- log(mammals$body)
Yfull <- log(mammals$brain)
nfull <- length(Yfull)

## test Gibbs sampler
rinit <- function(){
  return(list(chain_state =c(rnorm(2), rexp(1)), current_target = 0))
}

get_kernels <- function(X, Y){
  n <- length(Y)
  mle <- as.numeric((solve(t(X) %*% X, t(X) %*% Y))[,1])
  inv_XtX <- solve(t(X) %*% X)

  rinit <- function(){
    return(list(chain_state =c(rnorm(2), rexp(1)), current_target = 0))
  }
  single_kernel <- function(x, current_target){
    # sample from multivariate Gaussian beta | sigma^2
    cov_ <- x[3] * inv_XtX
    x12_ <- fast_rmvnorm(1, mean = mle, covariance = cov_)
    x[1] <- x12_[1]
    x[2] <- x12_[2]
    # sample from inverse Gamma sigma^2 | beta
    residuals <- (Y - X %*% t(x12_))
    x[3] <- rigamma(1, n/2, sum((residuals)^2) / 2)
    return(list(chain_state = x, current_target = current_target))
  }
  coupled_kernel <- function(x1, x2, current_target1, current_target2){
    cov_1 <- x1[3] * inv_XtX
    cov_2 <- x2[3] * inv_XtX
    x12_ <- gaussian_max_coupling(mle, mle, cov_1, cov_2) # outputs cbind of columns
    x1[1] <- x12_[1,1]
    x1[2] <- x12_[2,1]
    x2[1] <- x12_[1,2]
    x2[2] <- x12_[2,2]
    # sample from inverse Gamma sigma^2 | beta
    residuals1 <- (Y - X %*% x12_[,1,drop=F])
    residuals2 <- (Y - X %*% x12_[,2,drop=F])
    x3_ <- rigamma_coupled(n/2,n/2, sum((residuals1)^2) / 2, sum((residuals2)^2) / 2)
    x1[3] <- x3_[1]
    x2[3] <- x3_[2]
    return(list(chain_state1 = x1, chain_state2 = x2, current_target1 = current_target1, current_target2 = current_target2, coupled_prop = FALSE))
  }
  return(list(rinit = rinit, single_kernel = single_kernel, coupled_kernel = coupled_kernel))
}

# set training set size
n_T <- nfull/2
# set validation set size
n_V <- nfull - n_T
#
get_cv_estimator <- function(k, m){
  # draw random ordering of {1,...,n}
  random_order <- sample(1:nfull, size = nfull, replace = FALSE)
  # split data into training and validation sets
  Y_T <- (Yfull[random_order])[1:n_T]
  X_T <- (Xfull[random_order])[1:n_T]
  Y_V <- (Yfull[random_order])[(n_T+1):nfull]
  X_V <- (Xfull[random_order])[(n_T+1):nfull]
  kernels_ <- get_kernels(cbind(1,X_T), Y_T)
  testfunction <- function(x)  n_V * x[3] + sum((x[1] + x[2]*X_V - Y_V)^2)
  res <- unbiased_estimator(kernels_$single_kernel, kernels_$coupled_kernel, kernels_$rinit, h = testfunction, k = k, m = m)
  return(list(meetingtime = res$meetingtime, estimator = res$uestimator))
}

## test
# get_cv_estimator(k = 15, m = 25)

nrep <- 1000
results <- foreach(irep = 1:nrep) %dorng% {
  get_cv_estimator(k = 10, m = 25)
}
# qplot(x = sapply(results, function(x) x$meetingtime), geom = "histogram")
cv.errors <- sapply(results, function(x) x$estimator[[1]])
# summary(cv.errors)
mean(cv.errors)
sd(cv.errors) / sqrt(nrep)
qplot(x = cv.errors, geom = "histogram")

mean(cv.errors) - 2 * sd(cv.errors) / sqrt(nrep)
mean(cv.errors) + 2 * sd(cv.errors) / sqrt(nrep)

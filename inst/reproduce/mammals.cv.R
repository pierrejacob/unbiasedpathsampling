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

get_kernels <- function(X_T, Y_T, X_V, Y_V, theta){
  n_T <- length(Y_T)
  n_V <- length(Y_V)
  XX_T <- t(X_T) %*% X_T
  XX_V <- t(X_V) %*% X_V
  XY_T <- t(X_T) %*% Y_T
  XY_V <- t(X_V) %*% Y_V

  rinit <- function(){
    return(list(chain_state =c(rnorm(2), rexp(1)), current_target = 0))
  }

  single_kernel <- function(x, current_target){
    # sample from multivariate Gaussian beta | sigma^2
    lambda_ <- (XX_T + theta * XX_V) / x[3]
    cov_ <- solve(lambda_)
    mu_ <- cov_ %*% ((XY_T + theta * XY_V) / x[3])
    x12_ <- fast_rmvnorm(1, mean = mu_, covariance = cov_)
    x[1] <- x12_[1]
    x[2] <- x12_[2]
    # sample from inverse Gamma sigma^2 | beta
    residuals_T <- (Y_T - X_T %*% t(x12_))
    residuals_V <- (Y_V - X_V %*% t(x12_))
    x[3] <- rigamma(1, (n_T + theta * n_V)/2, sum((residuals_T^2)) / 2 + theta * sum((residuals_V^2)) / 2)
    return(list(chain_state = x, current_target = current_target))
  }

  coupled_kernel <- function(x1, x2, current_target1, current_target2){
    lambda_1 <- (XX_T + theta * XX_V) / x1[3]
    cov_1 <- solve(lambda_1)
    mu_1 <- cov_1 %*% ((XY_T + theta * XY_V) / x1[3])
    lambda_2 <- (XX_T + theta * XX_V) / x2[3]
    cov_2 <- solve(lambda_2)
    mu_2 <- cov_2 %*% ((XY_T + theta * XY_V) / x2[3])
    x12_ <- gaussian_max_coupling(mu_1, mu_2, cov_1, cov_2) # outputs cbind of columns
    x1[1] <- x12_[1,1]
    x1[2] <- x12_[2,1]
    x2[1] <- x12_[1,2]
    x2[2] <- x12_[2,2]
    # sample from inverse Gamma sigma^2 | beta
    residuals_T1 <- (Y_T - X_T %*% x12_[,1,drop=F])
    residuals_V1 <- (Y_V - X_V %*% x12_[,1,drop=F])
    residuals_T2 <- (Y_T - X_T %*% x12_[,2,drop=F])
    residuals_V2 <- (Y_V - X_V %*% x12_[,2,drop=F])
    x3_ <- rigamma_coupled((n_T + theta * n_V)/2, (n_T + theta * n_V)/2,
                           sum((residuals_T1^2)) / 2 + theta * sum((residuals_V1^2)) / 2,
                           sum((residuals_T2^2)) / 2 + theta * sum((residuals_V2^2)) / 2)
    x1[3] <- x3_[1]
    x2[3] <- x3_[2]
    return(list(chain_state1 = x1, chain_state2 = x2, current_target1 = current_target1, current_target2 = current_target2, coupled_prop = FALSE))
  }
  return(list(rinit = rinit, single_kernel = single_kernel, coupled_kernel = coupled_kernel))
}

n_T <- nfull/2
n_V <- nfull - n_T
get_cv_estimator <- function(k, m){
  # draw random ordering of {1,...,n}
  random_order <- sample(1:nfull, size = nfull, replace = FALSE)
  theta <- runif(1)
  # split data into training and validation sets
  Y_T <- (Yfull[random_order])[1:n_T]
  X_T <- (Xfull[random_order])[1:n_T]
  Y_V <- (Yfull[random_order])[(n_T+1):nfull]
  X_V <- (Xfull[random_order])[(n_T+1):nfull]
  kernels_ <- get_kernels(cbind(1,X_T), Y_T, cbind(1, X_V), Y_V, theta)
  testfunction <- function(x){
    return(-sum(dnorm(Y_V, mean = x[1] + x[2]*X_V, sd = sqrt(x[3]), log = TRUE)))
  }
  res <- unbiased_estimator(kernels_$single_kernel, kernels_$coupled_kernel, kernels_$rinit, h = testfunction, k = k, m = m)
  return(list(meetingtime = res$meetingtime, estimator = res$uestimator))
}
get_cv_estimator(k = 15, m = 25)

nrep <- 1000
results <- foreach(irep = 1:nrep) %dorng% {
  get_cv_estimator(k = 10, m = 25)
}
# summary(sapply(results, function(x) x$meetingtime))
# qplot(x = sapply(results, function(x) x$meetingtime), geom = "histogram")
cv.errors <- sapply(results, function(x) x$estimator[[1]])
# summary(cv.errors)
mean(cv.errors)
sd(cv.errors) / sqrt(nrep)
qplot(x = cv.errors, geom = "histogram")

mean(cv.errors) - 2 * sd(cv.errors) / sqrt(nrep)
mean(cv.errors) + 2 * sd(cv.errors) / sqrt(nrep)

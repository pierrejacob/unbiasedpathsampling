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

data("stackloss")

# head(stackloss)
# matplot(stackloss, type = "l")
stackloss.df <- stackloss
stackloss.df$index <- 1:nrow(stackloss.df)

g <- ggplot(reshape2::melt(stackloss.df, id = "index"), aes(x = index, y = value, group = variable, linetype = variable)) + geom_line()
g
# ggsave(filename = "stacklossdata.pdf", plot = g, width = 10, height = 5)

Y <- stackloss[,4]
X <- cbind(1, stackloss[,1:3])

Xfull <- as.matrix(X)
Yfull <- as.numeric(Y)
nfull <- length(Yfull)
ncovariates <- ncol(X)

## test Gibbs sampler

get_kernels <- function(X_T, Y_T, X_V, Y_V, theta){
  n_T <- length(Y_T)
  n_V <- length(Y_V)
  XX_T <- t(X_T) %*% X_T
  XX_V <- t(X_V) %*% X_V
  XY_T <- t(X_T) %*% Y_T
  XY_V <- t(X_V) %*% Y_V

  rinit <- function(){
    return(list(chain_state = c(rnorm(ncovariates), rexp(1)), current_target = 0))
  }

  single_kernel <- function(x, current_target){
    # sample from multivariate Gaussian beta | sigma^2
    lambda_ <- (XX_T + theta * XX_V) / x[ncovariates+1]
    cov_ <- solve(lambda_)
    mu_ <- cov_ %*% ((XY_T + theta * XY_V) / x[ncovariates+1])
    x12_ <- fast_rmvnorm(1, mean = mu_, covariance = cov_)
    for (icovariate in 1:ncovariates){
      x[icovariate] <-  x12_[icovariate]
    }
    # sample from inverse Gamma sigma^2 | beta
    residuals_T <- (Y_T - X_T %*% t(x12_))
    residuals_V <- (Y_V - X_V %*% t(x12_))
    x[ncovariates+1] <- rigamma(1, (n_T + theta * n_V)/2, sum((residuals_T^2)) / 2 + theta * sum((residuals_V^2)) / 2)
    return(list(chain_state = x, current_target = current_target))
  }

  coupled_kernel <- function(x1, x2, current_target1, current_target2){
    lambda_1 <- (XX_T + theta * XX_V) / x1[ncovariates+1]
    cov_1 <- solve(lambda_1)
    mu_1 <- cov_1 %*% ((XY_T + theta * XY_V) / x1[ncovariates+1])
    lambda_2 <- (XX_T + theta * XX_V) / x2[ncovariates+1]
    cov_2 <- solve(lambda_2)
    mu_2 <- cov_2 %*% ((XY_T + theta * XY_V) / x2[ncovariates+1])
    x12_ <- gaussian_max_coupling(mu_1, mu_2, cov_1, cov_2) # outputs cbind of columns
    for (icovariate in 1:ncovariates){
      x1[icovariate] <-  x12_[icovariate,1]
      x2[icovariate] <-  x12_[icovariate,2]
    }
    # sample from inverse Gamma sigma^2 | beta
    residuals_T1 <- (Y_T - X_T %*% x12_[,1,drop=F])
    residuals_V1 <- (Y_V - X_V %*% x12_[,1,drop=F])
    residuals_T2 <- (Y_T - X_T %*% x12_[,2,drop=F])
    residuals_V2 <- (Y_V - X_V %*% x12_[,2,drop=F])
    x3_ <- rigamma_coupled((n_T + theta * n_V)/2, (n_T + theta * n_V)/2,
                           sum((residuals_T1^2)) / 2 + theta * sum((residuals_V1^2)) / 2,
                           sum((residuals_T2^2)) / 2 + theta * sum((residuals_V2^2)) / 2)
    x1[ncovariates+1] <- x3_[1]
    x2[ncovariates+1] <- x3_[2]
    return(list(chain_state1 = x1, chain_state2 = x2, current_target1 = current_target1, current_target2 = current_target2, coupled_prop = FALSE))
  }
  return(list(rinit = rinit, single_kernel = single_kernel, coupled_kernel = coupled_kernel))
}

n_T <- nfull-1
n_V <- nfull - n_T

# split data into training and validation sets
# random_order <- sample(1:nfull, size = nfull, replace = FALSE)
# Y_T <- (Yfull[random_order])[1:n_T]
# Y_V <- (Yfull[random_order])[(n_T+1):nfull]
# X_T <- (Xfull[random_order,])[1:n_T,]
# X_V <- (Xfull[random_order,])[(n_T+1):nfull,,drop=F]
# theta <- 0.5
# kernels_ <- get_kernels(X_T, Y_T, X_V, Y_V, theta)
# state <- kernels_$rinit()
# niterations <- 1e3
# chain <- matrix(nrow = niterations, ncol = ncovariates+1)
# for (iteration in 1:niterations){
#   new_state <- kernels_$single_kernel(state$chain_state, state$current_target)
#   state <- new_state
#   chain[iteration,] <- state$chain_state
# }
# matplot(chain[,1], type = "l")
# matplot(chain[,2:4], type = "l")
# matplot(chain[,5], type = "l")
# colMeans(chain[100:niterations,1:4])
# lm(Yfull ~ Xfull[,2:4])

get_cv_estimator <- function(k, m){
  # draw random ordering of {1,...,n}
  random_order <- sample(1:nfull, size = nfull, replace = FALSE)
  theta <- runif(1)
  # split data into training and validation sets
  Y_T <- (Yfull[random_order])[1:n_T]
  Y_V <- (Yfull[random_order])[(n_T+1):nfull]
  X_T <- (Xfull[random_order,])[1:n_T,]
  X_V <- (Xfull[random_order,])[(n_T+1):nfull,,drop=F]
  kernels_ <- get_kernels(X_T, Y_T, X_V, Y_V, theta)
  testfunction <- function(x){
    return(-sum(dnorm(Y_V, mean = X_V %*% matrix(x[1:ncovariates], ncol = 1), sd = sqrt(x[ncovariates+1]), log = TRUE)))
  }
  res <- unbiased_estimator(kernels_$single_kernel, kernels_$coupled_kernel, kernels_$rinit, h = testfunction, k = k, m = m)
  return(list(meetingtime = res$meetingtime, estimator = res$uestimator[[1]], index = random_order[length(random_order)]))
}


nrep <- 10000
results <- foreach(irep = 1:nrep) %dorng% {
  get_cv_estimator(k = 10, m = 25)
}

indices <- sapply(results, function(x) x$index)
meetingtimes <- sapply(results, function(x) x$meetingtime)
# qplot(x = meetingtimes, geom = "blank") + geom_histogram(aes(y = ..density..))

estimators <- sapply(results, function(x) x$estimator)
# qplot(x = estimators, geom = "blank") + geom_histogram(aes(y = ..density..))

g <- qplot(x = indices, y = estimators, geom = "violin", group = indices)
g <- g + xlab("observation left out") + ylab("CV objective")
g
# ggsave(filename = "stackloss.cv.pdf", plot = g, width = 5, height = 5)

# qplot(x = indices, y = meetingtimes, geom = "violin", group = indices)

mean(estimators) - 2 * sd(estimators) / sqrt(nrep)
mean(estimators) + 2 * sd(estimators) / sqrt(nrep)



# Unbiased path sampling (UPS) between two multivariate Normal distributions.
# This is the multivariate version of inst/normal.R.
# pi_0 = Normal(0, I) and pi_1 = Normal(D * 1, I) in dimension 'dimension'.
# Both unnormalized densities have the same normalizing constant,
# so the target of UPS is log(Z_1 / Z_0) = log(1) = 0.

# load package
library(unbiasedpathsampling)
# load packages for parallel computation
library(doParallel)
library(doRNG)
# load package for data manipulation
library(dplyr)
# load packages for plotting
library(ggplot2)
# load graphical theme
setmytheme()
# register processing cores
registerDoParallel(cores = max(1, detectCores()-2))
# clean environment
rm(list = ls())
# set seed
set.seed(21)

## set timer
library(tictoc)
tic("UPS multivariate normal")


# dimension of the state space
dimension <- 5
# set path of distributions
D <- 4
mean1 <- rep(D, dimension)
# initial (log-density, up to a constant)
target0 <- function(x) -0.5 * sum(x^2)
# final (log-density, up to the same constant)
target1 <- function(x) -0.5 * sum((x - mean1)^2)
# intermediate
target_lambda <- function(x, lambda) (1 - lambda) * target0(x) + lambda * target1(x)
# derivative with respect to lambda
loggammaratio <- function(x) target1(x) - target0(x)
# For reference: pi_lambda = Normal(lambda * mean1, I) and
# E_lambda[loggammaratio(X)] = (lambda - 1/2) * sum(mean1^2), which integrates to zero over [0,1]
exact_E <- function(lambda) (lambda - 0.5) * sum(mean1^2)

# initial distribution for the MH chains
rinit <- function(lambda){
  chain_state <- rnorm(dimension, mean = -1, sd = 2)
  current_pdf <- target_lambda(chain_state, lambda)
  return(list(chain_state = chain_state, current_pdf = current_pdf))
}

# Random walk MH kernels (single and coupled) targeting pi_lambda,
# with Normal proposals of covariance Sigma_proposal (a vector here, i.e. diagonal covariance).
# The coupled kernel uses reflection-maximal couplings of the proposals, see ?get_mrth_kernels
get_kernels <- function(lambda, Sigma_proposal){
  kernels <- get_mrth_kernels(logtarget = function(x) target_lambda(x, lambda), Sigma_proposal = Sigma_proposal)
  return(list(sk = kernels$single_kernel, ck = kernels$coupled_kernel))
}

## Now let's find some reasonable values for k and m on a grid of values of lambda
lambda_grid <- 0:10/10
nrep <- 100
k_grid <- rep(0, length(lambda_grid))
# proposal variance: usual scaling 2.38^2/d for random walk MH, for a target with unit variances
Sigma_proposal <- rep(2.38^2 / dimension, dimension)
meetings.df <- data.frame()
for (ilambda in 1:length(lambda_grid)){
  lambda <- lambda_grid[ilambda]
  ri <- function() rinit(lambda)
  kernels <- get_kernels(lambda, Sigma_proposal)
  meetings <- foreach(i = 1:nrep, .combine = c) %dorng% {
    unbiased_estimator(kernels$sk, kernels$ck, ri)$meetingtime
  }
  k_grid[ilambda] <- as.numeric(floor(quantile(meetings, probs = 0.99)))
  meetings.df <- rbind(meetings.df, data.frame(ilambda = ilambda, lambda = lambda, meetings = meetings))
}
## plot meeting times for different lambdas
g <- ggplot(meetings.df, aes(x = lambda, y = meetings, group = ilambda)) + geom_violin() +
  geom_line(data=data.frame(x = lambda_grid, y = k_grid), aes(x = x, y = y, group = NULL))
g <- g + xlab(expression(lambda)) + ylab("meeting times")
g
# ggsave(filename = "ups.mvnormal.meetings.pdf", plot = g, width = 5, height = 5)

meantau <- as.numeric((meetings.df %>% group_by(ilambda, lambda) %>% summarise(meantau = mean(meetings)))$meantau)
# to have approximately constant cost
# we set m to be the maximal value of 5 * k
# and then we adjust by adding max E[tau] - E[tau] to each value
m_grid <- floor(5*max(k_grid) + max(meantau) - meantau)
# this is such that the expected cost m + E[tau] is roughly constant
m_grid + meantau

# Now, find first moments and second moments of each component
# of the distribution pi_lambda for lambda on the grid
m1s <- matrix(0, nrow = length(lambda_grid), ncol = dimension)
m2s <- matrix(0, nrow = length(lambda_grid), ncol = dimension)
averagecost <- rep(0, length(lambda_grid))
testfunctions <- list(h = function(x) x, h2 = function(x) x^2)
for (ilambda in 1:length(lambda_grid)){
  lambda <- lambda_grid[ilambda]
  ri <- function() rinit(lambda)
  kernels <- get_kernels(lambda, Sigma_proposal)
  uestimators_ <- foreach(i = 1:nrep) %dorng% {
    unbiased_estimator(kernels$sk, kernels$ck, ri, h = testfunctions, k = k_grid[ilambda], m = m_grid[ilambda])
  }
  m1_estimators <- t(sapply(uestimators_, function(x) x$uestimator[[1]]))
  m2_estimators <- t(sapply(uestimators_, function(x) x$uestimator[[2]]))
  cost_estimators <- sapply(uestimators_, function(x) 2*x$meetingtime + max(1, x$iteration - x$meetingtime + 1))
  m1s[ilambda,] <- colMeans(m1_estimators)
  m2s[ilambda,] <- colMeans(m2_estimators)
  averagecost[ilambda] <- mean(cost_estimators)
}
# compare estimated means with exact means lambda * D
cbind(lambda = lambda_grid, estimated_mean = rowMeans(m1s), exact_mean = lambda_grid * D)

# from these first and second moments, we can fit a better initial distribution
# and different tuning parameters for the MCMC kernel
# (here we use the estimated marginal variances, i.e. diagonal covariance matrices)
get_tuned_settings <- function(index_lambda){
  mean_target <- m1s[index_lambda,]
  # guard against negative variance estimates
  var_target <- pmax(m2s[index_lambda,] - m1s[index_lambda,]^2, 1e-2)
  return(list(mean_target = mean_target, var_target = var_target,
              Sigma_proposal = 2.38^2 / dimension * var_target,
              k = k_grid[index_lambda], m = m_grid[index_lambda]))
}

uestimator_givenlambda <- function(lambda){
  # find closest value in the grid
  index_lambda <- which.min(abs(lambda_grid - lambda))
  # use parameters associated with closest value
  settings <- get_tuned_settings(index_lambda)
  ri <- function(){
    chain_state <- rnorm(dimension, mean = settings$mean_target, sd = sqrt(settings$var_target))
    current_pdf <- target_lambda(chain_state, lambda)
    return(list(chain_state = chain_state, current_pdf = current_pdf))
  }
  kernels <- get_kernels(lambda, settings$Sigma_proposal)
  ue_ <- unbiased_estimator(kernels$sk, kernels$ck, ri, h = loggammaratio, k = settings$k, m = settings$m)
  return(ue_)
}

# Now that we have a grid, we can sample lambda from a prior q(lambda)
# and for each draw, look up the k and m obtained for the nearest lambda in the grid
# then, with that choice of k and m, we can get unbiased estimators of E_lambda[d/dlambda log gamma ratio(X)] for each lambda
UPS <- function(rp, dp, n){
  lambdas <- rp(n)
  results.df <- foreach(i = 1:n, .combine = rbind) %dorng% {
    lambda <- lambdas[i]
    unbiased_result <- uestimator_givenlambda(lambda)
    uestimator <- unbiased_result$uestimator[[1]] / dp(lambda) # unbiased estimator
    meetingtime <- unbiased_result$meetingtime # meeting time tau
    iteration <- unbiased_result$iteration # total number of iteration (max(m,tau))
    cost <- 2*meetingtime + max(1, iteration + 1 - meetingtime)
    data.frame(i = i, lambda = lambda, uestimator = uestimator, meetingtime = meetingtime, iteration = iteration, cost = cost)
  }
  return(results.df)
}

# What about optimal distributions for lambda?
# For each lambda in the grid, estimate the square root of the second moment of the
# unbiased estimator of E_lambda[d/dlambda log gamma ratio(X)]
nrep <- 100
esquare <- rep(0, length(lambda_grid))
for (ilambda in 1:length(lambda_grid)){
  lambda <- lambda_grid[ilambda]
  esquare_estimator <- foreach(i = 1:nrep, .combine = c) %dorng% {
    unbiased_result <- uestimator_givenlambda(lambda)
    unbiased_result$uestimator[[1]]^2
  }
  esquare[ilambda] <- sqrt(mean(esquare_estimator))
}

g <- ggplot(data.frame(lambda = lambda_grid, esquare = esquare), aes(x = lambda, y = esquare)) +
  geom_line() + geom_point()
g <- g + xlab(expression(lambda)) + ylab(expression(sqrt("second moment")))
g
# ggsave(filename = "ups.mvnormal.sqrtmeansquare.pdf", plot = g, width = 5, height = 5)

# piecewise constant density for lambda, proportional to sqrt(second moment)
nlambda_grid <- length(lambda_grid)
heights <- rep(0, nlambda_grid-1)
widths <- rep(0, nlambda_grid-1)
masses <- rep(0, nlambda_grid-1)
for (i in 1:(nlambda_grid-1)){
  heights[i] <- (esquare[i] + esquare[i+1]) / 2
  widths[i] <- lambda_grid[i+1] - lambda_grid[i]
  masses[i] <- widths[i] * heights[i]
}
heights <- heights / sum(masses)
masses <- masses / sum(masses)

r_ <- function(n){
  allocation <- sample(x = 1:(nlambda_grid-1), size = n, replace = TRUE, prob = masses)
  return(runif(n, min = (lambda_grid[allocation]), max = (lambda_grid[allocation] + widths[allocation])))
}
d_ <- function(x){
  index <- findInterval(x, lambda_grid, rightmost.closed = TRUE, all.inside = TRUE)
  return(heights[index])
}

## Test r_ and d_
# hist(r_(1e6), breaks = lambda_grid, prob = TRUE)
# curve(d_(x), add = TRUE, col = "red", n = 1000)

# Now path sampling with this proposal
noptimal <- 10000
results.optimal.df <- UPS(r_, d_, noptimal)
g <- ggplot(results.optimal.df, aes(x = lambda, y = uestimator)) + geom_point(alpha = 0.25)
g <- g + xlab(expression(lambda)) + ylab(expression(hat(E)(lambda)/q(lambda)))
g
# ggsave(filename = "ups.mvnormal.qopt.png", plot = g, width = 5, height = 5)

results.optimal.df %>% summarise(mean_cost = mean(cost), var_ups = var(uestimator))

## Come up with confidence interval at 95% for log(Z_1/Z_0), which is equal to 0 here
logZ1overZ0_hat_opt <- mean(results.optimal.df$uestimator)
logZ1overZ0_variance_opt <- var(results.optimal.df$uestimator)
logZ1overZ0_cost_opt <- mean(results.optimal.df$cost)
logZ1overZ0_error_opt <- sqrt(logZ1overZ0_variance_opt/nrow(results.optimal.df))
logZ1overZ0_CI_opt <- c(logZ1overZ0_hat_opt - 1.96*logZ1overZ0_error_opt, logZ1overZ0_hat_opt + 1.96*logZ1overZ0_error_opt)
cat("[", round(logZ1overZ0_CI_opt[1],2), ",", round(logZ1overZ0_CI_opt[2],2), "]\n")

## print time taken
toc()

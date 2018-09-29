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

# Set distributions
target0 <- function(x) - (x[1]+2)^2 - 0.5 * x[2]^2
target1 <- function(x) -0.1 * ( ((x[1]-1)^2-x[2]^2 )^2 + 10*(x[1]^2-5)^2+ (x[1]+x[2])^4 + (x[1]-x[2])^4)
# compute constant exactly
target0marginal1 <- function(x) integrate(function(a) sapply(a, function(z) exp(target0(c(x, z)))), lower = -5, upper = 5, subdivisions = 1e4)$value
Z0 <- integrate(function(x) sapply(x, function(z) target0marginal1(z)), lower = -5, upper = 5, subdivisions = 1e4)$value
target1marginal1 <- function(x) integrate(function(a) sapply(a, function(z) exp(target1(c(x, z)))), lower = -5, upper = 5, subdivisions = 1e4)$value
Z1 <- integrate(function(x) sapply(x, function(z) target1marginal1(z)), lower = -5, upper = 5, subdivisions = 1e4)$value
log(Z1/Z0)

## Plot initial distribution and target
xseq <- seq(from = -3, to = 3, length.out = 50)
yseq <- seq(from = -3, to = 3, length.out = 50)
df.target <- expand.grid(xseq, yseq)
df.target$target0 <- apply(df.target, MARGIN = 1, function(x) target0(x[1:2]))
df.target$target1 <- apply(df.target, MARGIN = 1, function(x) target1(x[1:2]))
ggplot(df.target, aes(x = Var1, y = Var2, z = target0)) + geom_contour(bins = 20)
ggplot(df.target, aes(x = Var1, y = Var2, z = target1)) + geom_contour(bins = 20)

# path of distributions
target_lambda <- function(x, lambda) (1 - lambda) * target0(x) + lambda * target1(x)
# derivative wrt lambda
loggammaratio <- function(x) target1(x) - target0(x)
# dimension of target
dimension <- 2
# initial distribution of chains
rinit <- function(lambda){
  chain_state <- fast_rmvnorm(1, mean = rep(-2, dimension), diag(1, dimension, dimension))[1,]
  current_pdf <- target_lambda(chain_state, lambda)
  return(list(chain_state = chain_state, current_pdf = current_pdf))
}

# Metropolis hastings for a certain lambda between 0 and 1
Sigma_q <- diag(2, dimension, dimension)
Sigma_chol <- chol(Sigma_q)
Sigma_chol_inv <- solve(Sigma_chol)

single_kernel <- function(chain_state, current_pdf, lambda){
  proposal_value <- chain_state + fast_rmvnorm_chol(1, rep(0, dimension), Sigma_chol)[1,]
  proposal_pdf <- target_lambda(proposal_value, lambda)
  accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
  if (accept){
    return(list(chain_state = proposal_value, current_pdf = proposal_pdf))
  } else {
    return(list(chain_state = chain_state, current_pdf = current_pdf))
  }
}


# Markov kernel of the coupled chain
coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, lambda){
  cproposals_ <- unbiasedpathsampling:::rnorm_reflectionmax_(chain_state1, chain_state2, Sigma_chol, Sigma_chol_inv)
  proposal1 <- cproposals_$xy[,1]
  proposal2 <- cproposals_$xy[,2]
  proposal_pdf1 <- target_lambda(proposal1, lambda)
  proposal_pdf2 <- target_lambda(proposal2, lambda)
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

lambda_grid <- 0:10/10
nrep <- 1000
k_grid <- rep(0, length(lambda_grid))
meetings.df <- data.frame()
for (ilambda in 1:length(lambda_grid)){
  lambda <- lambda_grid[ilambda]
  ri <- function() rinit(lambda)
  sk <- function(x, f) single_kernel(x, f, lambda)
  ck <- function(x1, x2, f1, f2) coupled_kernel(x1, x2, f1, f2, lambda)
  meetings <- foreach(i = 1:nrep, .combine = c) %dorng% {
    unbiased_estimator(sk, ck, ri)$meetingtime
  }
  k_grid[ilambda] <- 2*as.numeric(floor(quantile(meetings, probs = 0.99)))
  meetings.df <- rbind(meetings.df, data.frame(ilambda = ilambda, lambda = lambda, meetings = meetings))
}

g <- ggplot(meetings.df, aes(x = lambda, y = meetings, group = ilambda)) + geom_violin() +
  geom_line(data=data.frame(x = lambda_grid, y = k_grid), aes(x = x, y = y, group = NULL))
g <- g + xlab(expression(lambda)) + ylab("meeting times")
g
# ggsave(filename = "ups.doublewell.meetings.pdf", plot = g, width = 5, height = 5)

meantau <- as.numeric((meetings.df %>% group_by(ilambda, lambda) %>% summarise(meantau = mean(meetings)))$meantau)
# to have approximately constant cost
# we set m to be the maximal value of 5 * k
# and then we adjust by adding max E[tau] - E[tau] to each value
m_grid <- floor(5*max(k_grid) + max(meantau) - meantau)
# this is such that the expected cost m + E[tau] is roughly constant
m_grid + meantau

uestimator_givenlambda <- function(lambda){
  # index_lambda <- which(order(c(lambda_grid, lambda)) == (length(lambda_grid)+1))
  index_lambda <- which.min(abs(lambda_grid - lambda))
  k <- k_grid[index_lambda]
  m <- m_grid[index_lambda]
  ri <- function() rinit(lambda)
  sk <- function(x, f) single_kernel(x, f, lambda)
  ck <- function(x1, x2, f1, f2) coupled_kernel(x1, x2, f1, f2, lambda)
  ue_ <- unbiased_estimator(sk, ck, ri, h = loggammaratio, k = k, m = m)
  return(ue_)
}

# Now that we have a grid, we can sample lambda from a prior p(lambda)
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
    data.frame(i = 1, lambda = lambda, uestimator = uestimator, meetingtime = meetingtime, iteration = iteration, cost = cost)
  }
  return(results.df)
}

## Now fit a better proposal
nrep <- 100
esquare <- rep(0, length(lambda_grid))
for (ilambda in 1:length(lambda_grid)){
  lambda <- lambda_grid[ilambda]
  esquare_estimator <- foreach(i = 1:nrep, .combine = c) %dorng% {
    uestimator_givenlambda(lambda)$uestimator[[1]]^2
  }
  esquare[ilambda] <- sqrt(mean(esquare_estimator))
}

g <- qplot(x = lambda_grid, y = esquare, geom = "line") + geom_point() + ylim(0,32)
g <- g + xlab(TeX("$\\lambda$")) + ylab(TeX("$\\sqrt{second\\, moment}$"))
g

# ggsave(filename = "ups.doublewell.sqrtmeansquare.pdf", plot = g, width = 5, height = 5)

nlambda_grid <- length(lambda_grid)
heights <- rep(0, nlambda_grid-1)
widths <- rep(0, nlambda_grid-1)
masses <- rep(0, nlambda_grid-1)
for (i in 1:(nlambda_grid-1)){
  heights[i] <- (esquare[i] + esquare[i+1]) / 2
  widths[i] <- lambda_grid[i+1] - lambda_grid[i]
  masses[i] <- widths[i] * heights[i]
}
masses <- masses / sum(masses)

r_ <- function(n){
  allocation <- sample(x = 1:(nlambda_grid-1), size = n, replace = TRUE, prob = masses)
  return(runif(n, min = (lambda_grid[allocation]), max = (lambda_grid[allocation] + widths[allocation])))
}

d__ <- function(x){
  index <- which(order(c(x, lambda_grid)) == 1)
  if (index == 1){
    return(heights[1])
  } else {
    return(heights[index-1])
  }
}
d_integral <- integrate(f = function(x) sapply(x, d__), lower = 0, upper = 1, subdivisions = 1e3)
d_ <- function(x) sapply(x, d__)/d_integral$value

# Now path sampling with this proposal
noptimal <- 1000
results.optimal.df <- UPS(r_, d_, noptimal)
results.optimal.df %>% head

g <- ggplot(results.optimal.df, aes(x = lambda, y = uestimator)) + geom_point(alpha = 0.25)
g <- g + xlab(expression(lambda)) + ylab(TeX("$\\hat{E}(\\lambda)/q(\\lambda)$"))
g
# ggsave(filename = "ups.doublewell.qopt.png", plot = g, width = 5, height = 5)

## Confidence interval at 95%
logZ1overZ0_hat_opt <- mean(results.optimal.df$uestimator)
logZ1overZ0_variance_opt <- var(results.optimal.df$uestimator)
logZ1overZ0_cost_opt <- mean(results.optimal.df$cost)
logZ1overZ0_error_opt <- sqrt(logZ1overZ0_variance_opt/nrow(results.optimal.df))
cat("[", logZ1overZ0_hat_opt-1.96*logZ1overZ0_error_opt, ",", logZ1overZ0_hat_opt+1.96*logZ1overZ0_error_opt, "]\n")
log(Z1/Z0)



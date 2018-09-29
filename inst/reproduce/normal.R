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
# set path of distributions
D <- 4
# initial
target0 <- function(x) -0.5 * x^2
# final
target1 <- function(x) -0.5 * (x-D)^2
# intermediate
target_lambda <- function(x, lambda) (1 - lambda) * target0(x) + lambda * target1(x)
# derivative with respect to lambda
loggammaratio <- function(x) target1(x) - target0(x)

# initial distribution for the MH chain
rinit <- function(lambda){
  chain_state <- rnorm(1, mean = -1, sd = 2)
  current_pdf <- target_lambda(chain_state, lambda)
  return(list(chain_state = chain_state, current_pdf = current_pdf))
}
# maximal coupling of two distributions p and q
# specified by rp,dp (generate from p, and evaluate log-density of p),
# likewise for q with rq,dq
max_coupling <- function(rp, dp, rq, dq){
  x <- rp()
  if (dp(x) + log(runif(1)) < dq(x)){
    return(c(x,x))
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- rq()
      reject <- (dq(y) + log(runif(1)) < dp(y))
    }
    return(c(x,y))
  }
}

# Metropolis hastings for a certain lambda between 0 and 1
single_kernel <- function(chain_state, current_pdf, lambda, sigmaq){
  proposal_value <- chain_state + sigmaq * rnorm(1)
  proposal_pdf <- target_lambda(proposal_value, lambda)
  accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
  if (accept){
    return(list(chain_state = proposal_value, current_pdf = proposal_pdf))
  } else {
    return(list(chain_state = chain_state, current_pdf = current_pdf))
  }
}

# Markov kernel of the coupled chain
coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, lambda, sigmaq){
  proposals <- max_coupling(rp = function() rnorm(1, chain_state1, sigmaq),
                            dp = function(x) dnorm(x, chain_state1, sigmaq, log = TRUE),
                            rq = function() rnorm(1, chain_state2, sigmaq),
                            dq = function(x) dnorm(x, chain_state2, sigmaq, log = TRUE))
  proposal1 <- proposals[1]
  proposal2 <- proposals[2]
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


## Now let's find some reasonable values for k and m on a grid of values of lambda
lambda_grid <- 0:10/10
nrep <- 100
k_grid <- rep(0, length(lambda_grid))
sigmaq <- 1
meetings.df <- data.frame()
for (ilambda in 1:length(lambda_grid)){
  lambda <- lambda_grid[ilambda]
  ri <- function() rinit(lambda)
  sk <- function(x, f) single_kernel(x, f, lambda, sigmaq)
  ck <- function(x1, x2, f1, f2) coupled_kernel(x1, x2, f1, f2, lambda, sigmaq)
  meetings <- foreach(i = 1:nrep, .combine = c) %dorng% {
    unbiased_estimator(sk, ck, ri)$meetingtime
  }
  k_grid[ilambda] <- as.numeric(floor(quantile(meetings, probs = 0.99)))
  meetings.df <- rbind(meetings.df, data.frame(ilambda = ilambda, lambda = lambda, meetings = meetings))
}
## plot meeting times for different lambdas
g <- ggplot(meetings.df, aes(x = lambda, y = meetings, group = ilambda)) + geom_violin() +
  geom_line(data=data.frame(x = lambda_grid, y = k_grid), aes(x = x, y = y, group = NULL))
g <- g + xlab(expression(lambda)) + ylab("meeting times")
g
# ggsave(filename = "ups.normal.meetings.pdf", plot = g, width = 5, height = 5)


meantau <- as.numeric((meetings.df %>% group_by(ilambda, lambda) %>% summarise(meantau = mean(meetings)))$meantau)
# to have approximately constant cost
# we set m to be the maximal value of 5 * k
# and then we adjust by adding max E[tau] - E[tau] to each value
m_grid <- floor(5*max(k_grid) + max(meantau) - meantau)
# this is such that the expected cost m + E[tau] is roughly constant
m_grid + meantau

# Now, find first and second moment of distribution pi_lambda for lambda on the grid
m1s <- rep(0, length(lambda_grid))
m2s <- rep(0, length(lambda_grid))
averagecost <- rep(0, length(lambda_grid))
testfunctions <- list(h = function(x) x, h2 = function(x) x^2)
sigmaq <- 1
for (ilambda in 1:length(lambda_grid)){
  lambda <- lambda_grid[ilambda]
  ri <- function() rinit(lambda)
  sk <- function(x, f) single_kernel(x, f, lambda, sigmaq)
  ck <- function(x1, x2, f1, f2) coupled_kernel(x1, x2, f1, f2, lambda, sigmaq)
  uestimators_ <- foreach(i = 1:nrep) %dorng% {
    unbiased_estimator(sk, ck, ri, h = testfunctions, k = k_grid[ilambda], m = m_grid[ilambda])
  }
  m1_estimators <- sapply(uestimators_, function(x) x$uestimator[[1]])
  m2_estimators <- sapply(uestimators_, function(x) x$uestimator[[2]])
  cost_estimators <- sapply(uestimators_, function(x) 2*x$meetingtime + max(1, x$iteration - x$meetingtime + 1))
  m1s[ilambda] <- mean(m1_estimators)
  m2s[ilambda] <- mean(m2_estimators)
  averagecost[ilambda] <- mean(cost_estimators)
}

# from these first and second moments, we can fit a better initial distribution
# and different tuning parameters for the MCMC kernel

uestimator_givenlambda <- function(lambda){
  # find closest value in the grid
  index_lambda <- which.min(abs(lambda_grid - lambda))
  # use parameters associated with closest value
  mean_target <- m1s[index_lambda]
  var_target <- m2s[index_lambda] - m1s[index_lambda]^2
  sd_target <- sqrt(var_target)
  k <- k_grid[index_lambda]
  m <- m_grid[index_lambda]
  ri <- function(){
    chain_state <- rnorm(1, mean = mean_target, sd = sd_target)
    current_pdf <- target_lambda(chain_state, lambda)
    return(list(chain_state = chain_state, current_pdf = current_pdf))
  }
  sk <- function(x, f) single_kernel(x, f, lambda, sd_target)
  ck <- function(x1, x2, f1, f2) coupled_kernel(x1, x2, f1, f2, lambda, sd_target)
  ue_ <- unbiased_estimator(sk, ck, ri, h = loggammaratio, k = k, m = m)
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
    data.frame(i = 1, lambda = lambda, uestimator = uestimator, meetingtime = meetingtime, iteration = iteration, cost = cost)
  }
  return(results.df)
}

# What about optimal distributions for lambda?
nrep <- 100
esquare <- rep(0, length(lambda_grid))
for (ilambda in 1:length(lambda_grid)){
  lambda <- lambda_grid[ilambda]
  mean_target <- m1s[ilambda]
  var_target <- m2s[ilambda] - m1s[ilambda]^2
  sd_target <- sqrt(var_target)
  k <- k_grid[ilambda]
  m <- m_grid[ilambda]
  ri <- function(){
    chain_state <- rnorm(1, mean = mean_target, sd = sd_target)
    current_pdf <- target_lambda(chain_state, lambda)
    return(list(chain_state = chain_state, current_pdf = current_pdf))
  }
  sk <- function(x, f) single_kernel(x, f, lambda, sd_target)
  ck <- function(x1, x2, f1, f2) coupled_kernel(x1, x2, f1, f2, lambda, sd_target)
  esquare_estimator <- foreach(i = 1:nrep, .combine = c) %dorng% {
    unbiased_result <- unbiased_estimator(sk, ck, ri, h = loggammaratio, k = k, m = m)
    unbiased_result$uestimator[[1]]^2
  }
  esquare[ilambda] <- sqrt(mean(esquare_estimator))
}

g <- qplot(x = lambda_grid, y = esquare, geom = "line") + geom_point()
g <- g + xlab(TeX("$\\lambda$")) + ylab(TeX("$\\sqrt{second\\, moment}$")) + ylim(0, 10)
g
# ggsave(filename = "ups.normal.sqrtmeansquare.pdf", plot = g, width = 5, height = 5)

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

## Test r_ and d_
# hist(r_(1e6), breaks = lambda_grid, prob = TRUE)
# curve(d_(x), add = TRUE, col = "red", n = 1000)

# Now path sampling with this proposal
noptimal <- 5000
results.optimal.df <- UPS(r_, d_, noptimal)
g <- ggplot(results.optimal.df, aes(x = lambda, y = uestimator)) + geom_point(alpha = 0.25)
g <- g + xlab(expression(lambda)) + ylab(TeX("$\\hat{E}(\\lambda)/q(\\lambda)$"))
g
# ggsave(filename = "ups.normal.qopt.png", plot = g, width = 5, height = 5)

## Come up with confidence interval at 95%
logZ1overZ0_hat_opt <- mean(results.optimal.df$uestimator)
logZ1overZ0_variance_opt <- var(results.optimal.df$uestimator)
logZ1overZ0_cost_opt <- mean(results.optimal.df$cost)
logZ1overZ0_error_opt <- sqrt(logZ1overZ0_variance_opt/nrow(results.optimal.df))
cat("[", logZ1overZ0_hat_opt-1.96*logZ1overZ0_error_opt, ",", logZ1overZ0_hat_opt+1.96*logZ1overZ0_error_opt, "]\n")

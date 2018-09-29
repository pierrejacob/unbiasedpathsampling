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

unbiasedestimator <- function(single_kernel, coupled_kernel, rinit, h = function(x) x, k = 0, m = 1, max_iterations = Inf){
  chain_state1 <- rinit()
  chain_state2 <- rinit()
  # mcmcestimator computes the sum of h(X_t) for t=k,...,m
  mcmcestimator <- h(chain_state1)
  dimh <- length(mcmcestimator)
  if (k > 0){
    mcmcestimator <- rep(0, dimh)
  }
  # correction computes the sum of min(1, (t - k + 1) / (m - k + 1)) * (h(X_{t+1}) - h(X_t)) for t=k,...,max(m, tau - 1)
  correction <- rep(0, dimh)
  sres1 <- single_kernel(chain_state1)
  chain_state1 <- sres1
  if (k == 0){
    correction <- correction + (min(1, (0 - k + 1)/(m - k + 1))) * (h(chain_state1) - h(chain_state2))
  }
  if (k <= 1 && m >= 1){
    mcmcestimator <- mcmcestimator + h(chain_state1)
  }
  iter <- 1
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    if (meet){
      sres1 <- single_kernel(chain_state1)
      chain_state1 <- sres1
      chain_state2 <- chain_state1
      if (k <= iter && iter <= m){
        mcmcestimator <- mcmcestimator + h(chain_state1)
      }
    } else {
      res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2)
      chain_state1 <- res_coupled_kernel$chain_state1
      chain_state2 <- res_coupled_kernel$chain_state2
      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
      if (k <= iter){
        if (iter <= m){
          mcmcestimator <- mcmcestimator + h(chain_state1)
        }
        correction <- correction + (min(1, (iter-1 - k + 1)/(m - k + 1))) * (h(chain_state1) - h(chain_state2))
      }
    }
    # stop after max(m, tau) steps
    if (iter >= max(meetingtime, m)){
      finished <- TRUE
    }
  }
  mcmcestimator <- mcmcestimator / (m - k + 1)
  uestimator <- mcmcestimator + correction
  return(list(mcmcestimator = mcmcestimator, correction = correction, uestimator = uestimator,
              meetingtime = meetingtime, iteration = iter, finished = finished))
}

## Create data
p <- 7 # number of covariates
theta_star <- (0:6)/10 # data-generating parameter
expit <- function(z) 1 / (1 + exp(-z)) # expit function
n <- 1000 # number of observations
# generate X and Y from the model
X <- cbind(rep(1, n), scale(matrix(rnorm(n*(p-1)), ncol = p-1)))
logitprobs <- apply(X = X, MARGIN = 1, FUN = function(row) sum(row * theta_star))
Y <- rbinom(n, 1, expit(logitprobs))
# prior
b <- matrix(0, nrow = p, ncol = 1)
B <- diag(10, p, p)

# the lambda parameters multiplies the covariate matrix X
logistic_precomputation <- function(Y, X, b, B, lambda = 1){
  invB <- solve(B)
  invBtimesb <- invB %*% b
  Ykappa <- matrix(Y - rep(0.5, length(Y)), ncol=1)
  XTkappa <- t(lambda * X) %*% Ykappa
  KTkappaplusinvBtimesb <- XTkappa + invBtimesb
  return(list(n=nrow(X), p=ncol(X), X=lambda*X, Y=Y, b=b, B=B,
              invB=invB, invBtimesb=invBtimesb, KTkappaplusinvBtimesb=KTkappaplusinvBtimesb))
}

# MCMC transition kernel
single_kernel <- function(chain_state, logistic_setting){
  zs <- abs(unbiasedpathsampling:::xbeta_(logistic_setting$X, t(chain_state)))
  w <- BayesLogit::rpg(logistic_setting$n, h=1, z=zs)
  res <- unbiasedpathsampling:::m_sigma_function_(w, logistic_setting$X, logistic_setting$invB, logistic_setting$KTkappaplusinvBtimesb)
  chain_state <- t(fast_rmvnorm_chol(1, res$m, res$Cholesky))
  return(chain_state)
}

# Coupled MCMC transition kernel
coupled_kernel <- function(chain_state1, chain_state2, logistic_setting, return_ws=FALSE){
  ws <- sample_w(chain_state1, chain_state2, logistic_setting$X)
  betas <- sample_beta(ws$w1, ws$w2, logistic_setting)
  if (all(ws$w1 == ws$w2)){
    betas$beta2 <- betas$beta1
  }
  if(!return_ws){
    return(list(chain_state1=cbind(betas$beta1), chain_state2=cbind(betas$beta2)))
  } else {
    return(list(chain_state1=cbind(betas$beta1), chain_state2=cbind(betas$beta2),w1=ws$w1,w2=ws$w2))
  }
}
# initial distribution of the chains
rinit <- function(){
  t(fast_rmvnorm(1, mean = b, covariance = B))
}

lambda <- 0.5
logistic_setting <- logistic_precomputation(Y, X, b, B, lambda = lambda)
sk <- function(x) single_kernel(x, logistic_setting)
ck <- function(x1, x2) coupled_kernel(x1, x2, logistic_setting)
unbiasedestimator(sk, ck, rinit, h = function(x) x, k = 0, m = 0, max_iterations = 1e3)$meetingtime

testfunction <- function(beta, lambda){
  xbetas <- unbiasedpathsampling:::xbeta_(X, beta)
  return(sum(xbetas * Y - xbetas * exp(lambda * xbetas) / (1 + exp(lambda * xbetas))))
}


## meeting times
lambda_grid <- 0:10/10
nrep <- 1000
k_grid <- rep(0, length(lambda_grid))
meetings.df <- data.frame()
for (ilambda in 1:length(lambda_grid)){
  lambda <- lambda_grid[ilambda]
  logistic_setting <- logistic_precomputation(Y, X, b, B, lambda)
  sk <- function(x) single_kernel(x, logistic_setting)
  ck <- function(x1, x2) coupled_kernel(x1, x2, logistic_setting)
  meetings <- foreach(i = 1:nrep, .combine = c) %dorng% {
    unbiasedestimator(sk, ck, rinit)$meetingtime
  }
  k_grid[ilambda] <- as.numeric(floor(quantile(meetings, probs = 0.99)))
  meetings.df <- rbind(meetings.df, data.frame(ilambda = ilambda, lambda = lambda, meetings = meetings))
}

g <- ggplot(meetings.df, aes(x = lambda, y = meetings, group = ilambda)) + geom_violin() +
  geom_line(data=data.frame(x = lambda_grid, y = k_grid), aes(x = x, y = y, group = NULL))
g <- g + xlab(expression(lambda)) + ylab("meeting times")
g <- g + ylim(0,15)
g
# ggsave(filename = "ups.pgg.meetings.pdf", plot = g, width = 5, height = 5)

meantau <- as.numeric((meetings.df %>% group_by(ilambda, lambda) %>% summarise(meantau = mean(meetings)))$meantau)
# to have approximately constant cost
# we set m to be the maximal value of 5 * k
# and then we adjust by adding max E[tau] - E[tau] to each value
m_grid <- floor(5*max(k_grid) + max(meantau) - meantau)
# this is such that the expected cost m + E[tau] is roughly constant
m_grid + meantau

uestimator_givenlambda <- function(lambda){
  index_lambda <- which.min(abs(lambda_grid - lambda))
  k <- k_grid[index_lambda]
  m <- m_grid[index_lambda]
  logistic_setting <- logistic_precomputation(Y, X, b, B, lambda)
  sk <- function(x) single_kernel(x, logistic_setting)
  ck <- function(x1, x2) coupled_kernel(x1, x2, logistic_setting)
  ue_ <- unbiasedestimator(sk, ck, rinit, h = function(x) testfunction(x, lambda), k = k, m = m)
  return(ue_)
}
#
# Now that we have a grid, we can sample lambda from a prior p(lambda)
# and for each draw, look up the k and m obtained for the nearest lambda in the grid
# then, with that choice of k and m, we can get unbiased estimators of E_lambda[d/dlambda log gamma ratio(X)] for each lambda
UPS <- function(rp, dp, n){
  lambdas <- rp(n)
  results.df <- foreach(i = 1:n, .combine = rbind) %dorng% {
    lambda <- lambdas[i]
    unbiased_result <- uestimator_givenlambda(lambda)
    uestimator <- unbiased_result$uestimator / dp(lambda) # unbiased estimator
    meetingtime <- unbiased_result$meetingtime # meeting time tau
    iteration <- unbiased_result$iteration # total number of iteration (max(m,tau))
    cost <- 2*meetingtime + max(1, iteration + 1 - meetingtime)
    data.frame(i = 1, lambda = lambda, uestimator = uestimator, meetingtime = meetingtime, iteration = iteration, cost = cost)
  }
  return(results.df)
}

nrep <- 100
esquare <- rep(0, length(lambda_grid))
for (ilambda in 1:length(lambda_grid)){
  lambda <- lambda_grid[ilambda]
  esquare_estimator <- foreach(i = 1:nrep, .combine = c) %dorng% {
    uestimator_givenlambda(lambda)$uestimator^2
  }
  esquare[ilambda] <- sqrt(mean(esquare_estimator))
}

g <- qplot(x = lambda_grid, y = esquare, geom = "line") + geom_point()
g <- g + xlab(TeX("$\\lambda$")) + ylab(TeX("$\\sqrt{second\\, moment}$"))
g
# ggsave(filename = "ups.pgg.sqrtmeansquare.pdf", plot = g, width = 5, height = 5)

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
# ggsave(filename = "ups.pgg.qopt.png", plot = g, width = 5, height = 5)

## Confidence interval at 95%
logZ1overZ0_hat_opt <- mean(results.optimal.df$uestimator)
logZ1overZ0_variance_opt <- var(results.optimal.df$uestimator)
logZ1overZ0_cost_opt <- mean(results.optimal.df$cost)
logZ1overZ0_error_opt <- sqrt(logZ1overZ0_variance_opt/nrow(results.optimal.df))
cat("[", logZ1overZ0_hat_opt-1.96*logZ1overZ0_error_opt, ",", logZ1overZ0_hat_opt+1.96*logZ1overZ0_error_opt, "]\n")

## Now with a new grid
lambda_grid <- exp(seq(from = -10, to = 0, length.out = 10))
nrep <- 1000
k_grid <- rep(0, length(lambda_grid))
meetings.df <- data.frame()
for (ilambda in 1:length(lambda_grid)){
  lambda <- lambda_grid[ilambda]
  logistic_setting <- logistic_precomputation(Y, X, b, B, lambda)
  sk <- function(x) single_kernel(x, logistic_setting)
  ck <- function(x1, x2) coupled_kernel(x1, x2, logistic_setting)
  meetings <- foreach(i = 1:nrep, .combine = c) %dorng% {
    unbiasedestimator(sk, ck, rinit)$meetingtime
  }
  k_grid[ilambda] <- as.numeric(floor(quantile(meetings, probs = 0.99)))
  meetings.df <- rbind(meetings.df, data.frame(ilambda = ilambda, lambda = lambda, meetings = meetings))
}

g <- ggplot(meetings.df, aes(x = lambda, y = meetings, group = ilambda)) + geom_violin() +
  geom_line(data=data.frame(x = lambda_grid, y = k_grid), aes(x = x, y = y, group = NULL))
g <- g + xlab(expression(lambda)) + ylab("meeting times")
g <- g + scale_x_log10() + ylim(0,15)
g
# ggsave(filename = "ups.pgg.meetings.newgrid.pdf", plot = g, width = 5, height = 5)

meantau <- as.numeric((meetings.df %>% group_by(ilambda, lambda) %>% summarise(meantau = mean(meetings)))$meantau)
# to have approximately constant cost
# we set m to be the maximal value of 5 * k
# and then we adjust by adding max E[tau] - E[tau] to each value
m_grid <- floor(5*max(k_grid) + max(meantau) - meantau)
# this is such that the expected cost m + E[tau] is roughly constant
m_grid + meantau

nrep <- 100
esquare <- rep(0, length(lambda_grid))
for (ilambda in 1:length(lambda_grid)){
  lambda <- lambda_grid[ilambda]
  esquare_estimator <- foreach(i = 1:nrep, .combine = c) %dorng% {
    uestimator_givenlambda(lambda)$uestimator^2
  }
  esquare[ilambda] <- sqrt(mean(esquare_estimator))
}

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
g <- g + scale_x_log10()
g
# ggsave(filename = "ups.pgg.qopt.newgrid.png", plot = g, width = 5, height = 5)

## Confidence interval
logZ1overZ0_hat_opt <- mean(results.optimal.df$uestimator)
logZ1overZ0_variance_opt <- var(results.optimal.df$uestimator)
logZ1overZ0_cost_opt <- mean(results.optimal.df$cost)
logZ1overZ0_error_opt <- sqrt(logZ1overZ0_variance_opt/nrow(results.optimal.df))
cat("[", logZ1overZ0_hat_opt-1.96*logZ1overZ0_error_opt, ",", logZ1overZ0_hat_opt+1.96*logZ1overZ0_error_opt, "]\n")



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

# loglikelihood_cpp(Y, X, theta_mle)
# xthetamle <- (X %*% t(t(theta_mle)))
# sum(Y * xthetamle - log(1+exp(xthetamle)))
# sum(xthetamle * Y - log(1+exp(xthetamle)))

library(numDeriv)
hess <- hessian(func = function(theta) - loglikelihood_cpp(Y, X, theta), x = theta_mle)
asymptotic_var <- solve(hess)
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




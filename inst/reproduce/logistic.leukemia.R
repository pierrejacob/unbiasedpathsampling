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
  # current_nsamples1 <- current_nsamples1 + 1
  # samples1[current_nsamples1,] <- chain_state1
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

## Load leukemia data
library(BGPhazard)
data("leukemiaFZ")
leukemiaFZ$success <- leukemiaFZ$time >= 50
Y <- leukemiaFZ$success
X <- cbind(leukemiaFZ$wbc, leukemiaFZ$AG)

p <- 2 # number of covariates
expit <- function(z) 1 / (1 + exp(-z)) # expit function
n <- nrow(X) # number of observations
b <- matrix(0, nrow = p, ncol = 1)
B <- diag(10, p, p)

# the theta parameters multiplies the covariate matrix X at rows given by indices
logistic_precomputation <- function(Y, X, b, B, indices = NA, theta = 1){
  Xtilde <- X
  if (!is.na(indices)){
    Xtilde[indices,] <- theta * Xtilde[indices,]
  }
  invB <- solve(B)
  invBtimesb <- invB %*% b
  Ykappa <- matrix(Y - rep(0.5, length(Y)), ncol=1)
  XTkappa <- t(Xtilde) %*% Ykappa
  KTkappaplusinvBtimesb <- XTkappa + invBtimesb
  return(list(n=nrow(X), p=ncol(X), X=Xtilde, Y=Y, b=b, B=B,
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

rinit <- function(){
  t(fast_rmvnorm(1, mean = b, covariance = B))
}

indices <- 1
theta <- 0.5
logistic_setting <- logistic_precomputation(Y, X, b, B, indices = indices, theta = theta)
sk <- function(x) single_kernel(x, logistic_setting)
ck <- function(x1, x2) coupled_kernel(x1, x2, logistic_setting)
unbiasedestimator(sk, ck, rinit, h = function(x) x, k = 0, m = 0, max_iterations = 1e3)$meetingtime

testfunction <- function(beta, indices, theta){
  xbetas <- unbiasedpathsampling:::xbeta_(X[indices,,drop=F], beta)
  return(sum(xbetas * Y[indices] - xbetas * exp(theta * xbetas) / (1 + exp(theta * xbetas))))
}


#
uestimator_givenindices <- function(indices, k = 0, m = 0){
  theta <- runif(1)
  logistic_setting <- logistic_precomputation(Y, X, b, B, indices, theta)
  sk <- function(x) single_kernel(x, logistic_setting)
  ck <- function(x1, x2) coupled_kernel(x1, x2, logistic_setting)
  ue_ <- unbiasedestimator(sk, ck, rinit, h = function(x) testfunction(x, indices, theta), k = k, m = m)
  ue_$theta <- theta
  return(ue_)
}

# get meeting times
nrep <- 1000
meetings <- foreach(irep = 1:nrep) %dorng% {
  indices <- sample(1:nrow(X), size = 1)
  ue <- uestimator_givenindices(indices)
  ue$indices <- indices
  ue
}

meetingtimes <- sapply(meetings, function(x) x$meetingtime)
# hist(meetingtimes)
indices <-  sapply(meetings, function(x) x$indices)
g <- qplot(x = indices, y = meetingtimes, group = indices, geom = "violin")
g <- g + xlab("observation left out") + ylab("meeting times")
g
# ggsave(filename = "ups.leukemia.meetingtimes.pdf", plot = g, width = 5, height = 5)

k <- 100
m <- 5*k

nrep <- 10000
leave1out <- foreach(irep = 1:nrep) %dorng% {
  indices <- sample(1:nrow(X), size = 1)
  ue <- uestimator_givenindices(indices, k, m)
  ue$indices <- indices
  ue
}

cv <- sapply(leave1out, function(x) x$uestimator)-log(2)
meetings <- sapply(leave1out, function(x) x$meetingtime)
summary(meetings)
indices <- sapply(leave1out, function(x) x$indices)
# summary(cv)
# hist(cv)
# sd(cv)/sqrt(nrep)
cat("[", mean(cv)-2*sd(cv)/sqrt(nrep), mean(cv)+2*sd(cv)/sqrt(nrep), "]\n")

g <- qplot(x = cv, geom = "blank") + geom_histogram(aes(y = ..density..))
g <- g + xlab("CV objective")
g
# ggsave(filename = "ups.leukemia.cv.pdf", plot = g, width = 5, height = 5)

g <- qplot(x = indices, y = cv, group = indices, geom = "violin")
g <- g + xlab("observation left out") + ylab("CV objective")
g
# ggsave(filename = "ups.leukemia.cvindex.pdf", plot = g, width = 5, height = 5)


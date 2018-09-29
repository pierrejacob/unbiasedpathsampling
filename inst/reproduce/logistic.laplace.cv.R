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
dimension <- p
theta_star <- (0:6)/10 # data-generating parameter
expit <- function(z) 1 / (1 + exp(-z)) # expit function
n <- 1000 # number of observations
# generate X and Y from the model
X <- cbind(rep(1, n), scale(matrix(rnorm(n*(p-1)), ncol = p-1)))
logitprobs <- apply(X = X, MARGIN = 1, FUN = function(row) sum(row * theta_star))
Y <- rbinom(n, 1, expit(logitprobs))
b <- matrix(0, nrow = p, ncol = 1)
B <- diag(10, p, p)

library(Rcpp)
cppFunction('
            double loglikelihood_cpp(const NumericVector & Y, const NumericMatrix & X, const NumericVector & beta){
            double logpr = 0.;
            int p = beta.size();
            int ny = X.rows();
            for (int i = 0; i < ny; i++){
            double betaX = 0;
            for (int j = 0; j < p; j++){
            betaX = betaX + beta(j) * X(i,j);
            }
            logpr = logpr - log(1+exp(-betaX)) + (1-Y(i)) * (-betaX);
            }
            return(logpr);
            }')

logisticmle <- optim(par = theta_star, fn = function(theta) - loglikelihood_cpp(Y, X, theta))
theta_mle <- logisticmle$par
library(numDeriv)
hess <- hessian(func = function(theta) - loglikelihood_cpp(Y, X, theta), x = theta_mle)
asymptotic_var <- solve(hess)

cppFunction('
            double loglikelihood_index_theta_cpp(const NumericVector & Y, const NumericMatrix & X, const NumericVector & beta,
int index, double theta){
            double logpr = 0.;
            int p = beta.size();
            int ny = X.rows();
            for (int i = 0; i < ny; i++){
              double betaX = 0;
              for (int j = 0; j < p; j++){
                betaX = betaX + beta(j) * X(i,j);
              }
              if (i == index){
                logpr = logpr + theta * (- log(1+exp(-betaX)) + (1-Y(i)) * (-betaX));
              } else {
                logpr = logpr - log(1+exp(-betaX)) + (1-Y(i)) * (-betaX);
              }
            }
            return(logpr);
            }')


## Now path from Normal approx to posterior
# target0 <- function(x) fast_dmvnorm(t(x), theta_mle, asymptotic_var)
# target1 <- function(x) loglikelihood_cpp(Y, X, x) + fast_dmvnorm(t(x), b, B)
index <- 1
target_theta <- function(x, theta) loglikelihood_index_theta_cpp(Y,X,x,index,theta) + fast_dmvnorm(t(x), b, B)
loggammaratio <- function(x, theta){
  return(sum(x * X[index,]) * Y[index] - log(1 + exp(sum(x * X[index,]))))
}

# initial distribution of the chains
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

uestimator_givenindices <- function(index){
  theta <- runif(1)
  ri <- function() rinit(theta)
  sk <- function(x, f) single_kernel(x, f, theta)
  ck <- function(x1, x2, f1, f2) coupled_kernel(x1, x2, f1, f2, theta)
  ue_ <- unbiased_estimator(sk, ck, ri, h = loggammaratio, k = k, m = m)
  return(ue_)
}


nrep <- 100
leave1out <- foreach(irep = 1:nrep) %dorng% {
  uestimator_givenindices(1)
}
leave1outestimators <- sapply(leave1out, function(x) x$uestimator[[1]])
# hist(as.numeric(leave1outestimators))
# summary(leave1outestimators)
# summary(sapply(leave1out, function(x) x$meetingtime))

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -0.2534 -0.2196 -0.2122 -0.2120 -0.2040 -0.1778

# now randomizing over the indices


nrep <- 1000
leave1out <- foreach(irep = 1:nrep) %dorng% {
  index <- sample(1:nrow(X), size = 1)
  target_theta <- function(x, theta) loglikelihood_index_theta_cpp(Y,X,x,index,theta) + fast_dmvnorm(t(x), b, B)
  loggammaratio <- function(x, theta){
    return(sum(x * X[index,]) * Y[index] - log(1 + exp(sum(x * X[index,]))))
  }
  uestimator_givenindices <- function(index){
    theta <- runif(1)
    ri <- function() rinit(theta)
    sk <- function(x, f) single_kernel(x, f, theta)
    ck <- function(x1, x2, f1, f2) coupled_kernel(x1, x2, f1, f2, theta)
    ue_ <- unbiased_estimator(sk, ck, ri, h = loggammaratio, k = k, m = m)
    return(ue_)
  }
  uestimator_givenindices(index)
}

cv <- sapply(leave1out, function(x) x$uestimator[[1]])
cat("[", mean(cv)-2*sd(cv)/sqrt(nrep), mean(cv)+2*sd(cv)/sqrt(nrep), "]\n")

g <- qplot(x = cv, geom = "blank") + geom_histogram(aes(y = ..density..))
g <- g + xlab("CV objective")
g
# ggsave(filename = "ups.logistic.otherpath.cv.pdf", plot = g, width = 5, height = 5)
#



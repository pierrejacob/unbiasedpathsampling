###
## This folder contains scripts to reproduce the figures of the technical report
## on unbiased estimation of logarithms of normalizing constants,
## by Maxime Rischard, Pierre E. Jacob, Natesh Pillai
###

## On top of the unbiasedpathsampling package
## the scripts use doParallel, doRNG, dplyr, ggplot2, and latex2exp (see below about BayesLogit).
## The package unbiasedpathsampling requires Rcpp and RcppEigen to compile.
## The data about leukemia is in BGPhazard. The stack loss data is in R by default.
## The data about mammal weights is in MASS.
## The Laplace approximation in the logistic regression example uses the package numDeriv.

## The compute times indicated below depend on the number of available cores.
## The reported numbers are obtained on a 2015 laptop using 6 cores at 2.2 GHz.

## The following description follows the order of the results in the paper.

## Section 3.1.1 Normal example (compute time about 1 minute)
# normal.R creates the three plots of Figure 1 and the confidence intervals
system.time(source("normal.R"))

## Section 3.1.2 Double-well example (compute time about 1 minute)
# doublewell.R creates the three plots of Figure 2 and the confidence intervals
system.time(source("doublewell.R"))

## Section 3.2 Logistic regression (total compute time about 10 minutes)
## The code of this section requires BayesLogit, available on the archives of CRAN:
## https://cran.r-project.org/src/contrib/Archive/BayesLogit/
## It can be easily installed with devtools:
# library(devtools)
# install_version("BayesLogit")

## Section 3.2.1 Normalizing constant estimation
# logistic.pgg.R creates the three plots of Figure 3, and the two first confidence intervals
system.time(source("logistic.pgg.R"))
# logistic.laplace.R gives the last confidence interval in Section 3.2.1,
# corresponding to a path from a Laplace approx to the posterior
system.time(source("logistic.laplace.R"))

## Section 3.2.2 Cross-validation
# logistic.pgg.cv.R creates the first confidence interval and Figure 4.a.
system.time(source("logistic.pgg.cv.R"))
# logistic.laplace.cv.R creates the second confidence interval and Figure 4.b.
system.time(source("logistic.laplace.cv.R"))

## Section 3.2.3 Leukemia survival data (requires the BGPhazard package)
# logistic.leukemia.R creates Figure 5 and the confidence interval
system.time(source("logistic.leukemia.R"))

## Section 3.3 Linear regressions
## Section 3.3.1 Mammal weight data (takes ~ 10 seconds)
# mammals.cv.mse.R creates the cross-validation results, when the point prediction is assessed via mean squared error
system.time(source("mammals.cv.mse.R"))
# mammals.cv.R creates the cross-validation results under the logarithmic scoring rule
system.time(source("mammals.cv.R"))

## Section 3.3.1 Stack loss data (takes ~ 1 minute)
# stackloss.cv.R creates Figure 6 and the confidence interval for cross-validation
system.time(source("stackloss.cv.R"))



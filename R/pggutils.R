
#'@export
sample_w <- function(beta1, beta2, X){
  w1w2_mat <- unbiasedpathsampling:::w_max_couplingC(beta1, beta2, X)
  w1w2 <- list(w1=w1w2_mat[,1], w2=w1w2_mat[,2])
  return(w1w2)
}


# from sample_beta --------------------------------------------------------
#'@export
sample_beta <- function(w1, w2, logistic_setting, mc_prob=0.5){
  KTkappaplusinvBtimesb <- logistic_setting$KTkappaplusinvBtimesb
  res1 <- unbiasedpathsampling:::m_sigma_function_(w1, logistic_setting$X, logistic_setting$invB, KTkappaplusinvBtimesb)
  res2 <- unbiasedpathsampling:::m_sigma_function_(w2, logistic_setting$X, logistic_setting$invB, KTkappaplusinvBtimesb)
  x <- unbiasedpathsampling:::gaussian_max_coupling_cholesky(res1$m, res2$m, res1$Cholesky, res2$Cholesky, res1$Cholesky_inverse, res2$Cholesky_inverse)
  beta1 <- x[,1]
  beta2 <- x[,2]
  return(list(beta1=beta1, beta2=beta2))
}

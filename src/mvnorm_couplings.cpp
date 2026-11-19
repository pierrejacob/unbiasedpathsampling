#include <RcppEigen.h>
#include "mvnorm.h"
using namespace Rcpp;

// Couplings of multivariate Normal distributions.
// Adapted from the 'unbiasedpoisson' package (src/mvnorm.cpp).
// Each function returns a list with two entries:
// "xy" is a matrix with two columns, one for x and one for y, which are Normally distributed vectors;
// "identical" is a boolean indicating whether x = y.
// The helpers rmvnorm, rmvnorm_cholesky, dmvnorm and dmvnorm_cholesky_inverse are in mvnorm.cpp.

// sample from maximal coupling of Normal(mu1, Sigma1) and Normal(mu2, Sigma2)
// [[Rcpp::export]]
Rcpp::List rmvnorm_max_coupling_(const NumericVector & mu1, const NumericVector & mu2,
                                 const NumericMatrix & Sigma1, const NumericMatrix & Sigma2){
  RNGScope scope;
  int p = Sigma1.cols();
  bool identical_ = false;
  NumericMatrix xy(p,2);
  NumericMatrix x1 = rmvnorm(1, mu1, Sigma1);
  NumericMatrix x2(1,p);
  std::fill(x2.begin(), x2.end(), 0.);
  int nattempts = 0;
  for (int index = 0; index < p; index++){
    xy(index,0) = x1(0,index);
  }
  double logu = log(runif(1)(0));
  if ((dmvnorm(x1, mu1, Sigma1)(0) + logu) < dmvnorm(x1, mu2, Sigma2)(0)){
    for (int index = 0; index < p; index++){
      xy(index,1) = xy(index,0);
    }
    identical_ = true;
  } else {
    bool accept = false;
    while (!accept){
      nattempts++;
      x2 = rmvnorm(1, mu2, Sigma2);
      if ((dmvnorm(x2, mu2, Sigma2)(0) + log(runif(1)(0))) > dmvnorm(x2, mu1, Sigma1)(0)){
        accept = true;
      }
    }
    for (int index = 0; index < p; index++){
      xy(index,1) = x2(0,index);
    }
  }
  return Rcpp::List::create(Named("xy") = xy, Named("identical") = identical_, Named("nattempts") = nattempts);
}

// sample from maximal coupling of Normal(mu1, Sigma1) and Normal(mu2, Sigma2),
// given upper-triangular Cholesky factors of Sigma1, Sigma2 (e.g. chol(Sigma) in R)
// and their inverses (e.g. solve(chol(Sigma)) in R)
// [[Rcpp::export]]
Rcpp::List rmvnorm_max_coupling_cholesky_(const NumericVector & mu1,
                                          const NumericVector & mu2,
                                          const Eigen::MatrixXd & Cholesky1,
                                          const Eigen::MatrixXd & Cholesky2,
                                          const Eigen::MatrixXd & Cholesky_inverse1,
                                          const Eigen::MatrixXd & Cholesky_inverse2){
  RNGScope scope;
  int p = Cholesky1.cols();
  NumericMatrix xy(p,2);
  NumericMatrix x2(1,p);
  std::fill(x2.begin(), x2.end(), 0.);
  bool identical_ = false;
  NumericMatrix x1 = rmvnorm_cholesky(1, mu1, Cholesky1);
  for (int index = 0; index < p; index++){
    xy(index,0) = x1(0,index);
  }
  double u1 = log(runif(1)(0));
  if ((dmvnorm_cholesky_inverse(x1, mu1, Cholesky_inverse1)(0) + u1) <= dmvnorm_cholesky_inverse(x1, mu2, Cholesky_inverse2)(0)){
    for (int index = 0; index < p; index++){
      xy(index,1) = xy(index,0);
    }
    identical_ = true;
  } else {
    bool accept = false;
    while (!accept){
      x2 = rmvnorm_cholesky(1, mu2, Cholesky2);
      double u2 = log(runif(1)(0));
      if ((dmvnorm_cholesky_inverse(x2, mu2, Cholesky_inverse2)(0) + u2) >
            dmvnorm_cholesky_inverse(x2, mu1, Cholesky_inverse1)(0)){
        accept = true;
      }
    }
    for (int index = 0; index < p; index++){
      xy(index,1) = x2(0,index);
    }
  }
  return Rcpp::List::create(Named("xy") = xy, Named("identical") = identical_);
}

// Draw from a reflection-maximal coupling of Normal(mu1, Sigma) and Normal(mu2, Sigma)
// As described in "Coupling and Convergence for Hamiltonian Monte Carlo" by Bou-Rabee et al, arXiv:1805.00452v1
// arguments Sigma_chol and inv_Sigma_chol can be obtained e.g. in R as
// chol(Sigma) and solve(chol(Sigma))
// [[Rcpp::export]]
Rcpp::List rmvnorm_reflection_max_coupling_(const Eigen::VectorXd & mu1, const Eigen::VectorXd & mu2,
                                            const Eigen::MatrixXd & Sigma_chol, const Eigen::MatrixXd & inv_Sigma_chol){
  RNGScope scope;
  int d = mu1.size();
  NumericMatrix xy(d,2);
  bool identical_ = false;
  Eigen::VectorXd scaled_diff = (mu2-mu1).transpose() * inv_Sigma_chol;
  // generate xi, a standard Normal draw
  Eigen::VectorXd xi = as<Eigen::VectorXd>(rnorm(d, 0., 1.));
  Eigen::VectorXd eta;
  // define z and e
  Eigen::VectorXd z = - scaled_diff;
  double normz = z.norm();
  if (normz < 1e-15){
    // mu1 and mu2 are (numerically) identical: return identical draws,
    // which avoids dividing by normz below
    identical_ = true;
    xi = mu1 + (xi.transpose() * Sigma_chol).transpose();
    for(int i=0; i<d; i++){
      xy(i,0) = xi(i);
      xy(i,1) = xi(i);
    }
  } else {
    Eigen::ArrayXd e = z.array() / normz;
    NumericVector utilde = runif(1);
    double edotxi = (e * xi.array()).sum();
    if (log(utilde(0)) < (-0.5 * (edotxi + normz) * (edotxi + normz) + 0.5 * edotxi * edotxi)){
      eta = (xi.array() + z.array()).matrix();
      identical_ = true;
    } else {
      eta = (xi.array() - 2. * edotxi * e).matrix();
    }
    // construct x ~ Normal(mu1, Sigma) and y ~ Normal(mu2, Sigma) from xi and eta
    xi = mu1 + (xi.transpose() * Sigma_chol).transpose();
    eta = mu2 + (eta.transpose() * Sigma_chol).transpose();
    for(int i=0; i<d; i++){
      xy(i,0) = xi(i);
      if (identical_){
        xy(i,1) = xi(i);
      } else {
        xy(i,1) = eta(i);
      }
    }
  }
  return Rcpp::List::create(Named("xy") = xy, Named("identical") = identical_);
}

// Same as above, in the case where the covariance matrix is diagonal;
// Sigma_chol is then the vector of standard deviations.
// [[Rcpp::export]]
Rcpp::List rmvnorm_reflection_max_coupling_diag_(const Eigen::ArrayXd & mu1, const Eigen::ArrayXd & mu2,
                                                 const Eigen::ArrayXd & Sigma_chol){
  RNGScope scope;
  int d = mu1.size();
  NumericMatrix xy(d,2);
  bool identical_ = false;
  Eigen::ArrayXd scaled_diff = (mu2-mu1) / Sigma_chol;
  // generate xi, a standard Normal draw
  Eigen::ArrayXd xi = as<Eigen::ArrayXd>(rnorm(d, 0., 1.));
  Eigen::ArrayXd eta;
  // define z and e
  Eigen::ArrayXd z = - scaled_diff;
  double normz = z.matrix().norm();
  if (normz < 1e-15){
    identical_ = true;
    xi = mu1 + xi * Sigma_chol;
    for(int i=0; i<d; i++){
      xy(i,0) = xi(i);
      xy(i,1) = xi(i);
    }
  } else {
    Eigen::ArrayXd e = z / normz;
    NumericVector utilde = runif(1);
    double edotxi = (e * xi).sum();
    if (log(utilde(0)) < (-0.5 * (edotxi + normz) * (edotxi + normz) + 0.5 * edotxi * edotxi)){
      eta = xi + z;
      identical_ = true;
    } else {
      eta = xi - 2. * edotxi * e;
    }
    // construct x ~ Normal(mu1, Sigma) and y ~ Normal(mu2, Sigma) from xi and eta
    xi = mu1 + xi * Sigma_chol;
    eta = mu2 + eta * Sigma_chol;
    for (int i=0; i<d; i++){
      xy(i,0) = xi(i);
      if (identical_){
        xy(i,1) = xi(i);
      } else {
        xy(i,1) = eta(i);
      }
    }
  }
  return Rcpp::List::create(Named("xy") = xy, Named("identical") = identical_);
}

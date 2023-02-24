// Contains functions for generating multivariate normal samples (using Eigen) and 
// Generalized Inverse Gaussian samples (using the algorithm of Hormann, J. (2014)).

#include <RcppEigen.h>
#include <R.h>
#include <cmath>
#include <algorithm>
#include "rand_functions.h"

double epsilon_gig = std::sqrt(std::numeric_limits<double>::min());

Eigen::MatrixXd rmvnorm(const int &n, const Eigen::MatrixXd &mean, const Eigen::MatrixXd &cov){
  // generate samples from a multivariate normal distribution
  int dim = cov.rows();
  Eigen::MatrixXd K(cov.llt().matrixL()); // Cholesky decomposition
  Eigen::MatrixXd samps = Eigen::MatrixXd::Zero(dim, n);
  
  for(int j=0; j<n; j++){
    for(int i=0; i<dim; i++){
      // Eigen is column major by default, so fill by column
      samps(i, j) = R::rnorm(0.0, 1.0);
    }
  }
  
  K*=samps;
  
  K.colwise() += mean.col(0);
  
  return K;
}

double rgig1(double const &lambda, double const &chi, double const &psi){
  // generate 1 sample from the Generalized Inverse Gaussian distribution
  // using the algorithm of Hormann, J. (2014).
  
  double lambda_abs = std::abs(lambda);
  double rv;
  double alpha = std::sqrt(psi/chi);
  double beta = std::sqrt(psi*chi);
  if(std::abs(chi) < epsilon_gig){
    // generate gamma RV
    if (lambda > 0.0){
      rv = R::rgamma(lambda_abs, 2/psi);
    } else {
      rv = 1.0/R::rgamma(lambda_abs, 2.0/psi);
    }
    return rv;
    
  } else if ((lambda_abs >= 0.0) && (lambda_abs < 1) && (beta > 0) && (beta <= 2*std::sqrt(1-lambda_abs)/3)){
    rv = algo1(lambda_abs, beta);
    
  } else if (lambda_abs > 1.0 || beta > 1.0) {
    rv = algo3(lambda_abs, beta);
    
  } else {
    rv = algo2(lambda_abs, beta);
  }
  
  if(lambda > 0){
    return rv/alpha;
  } else {
    return 1/(alpha*rv);
  }
  
}

double algo1(double const &lambda, double const &beta){
  double m = mode(lambda, beta);
  double x0 = beta/(1-lambda);
  double x_star = std::max(x0, 2.0/beta);
  double k1 = g(m, lambda, beta);
  double A1 = k1*x0;
  double A2, k2;
  if (x0 < 2.0/beta){
    k2 = std::exp(-beta);
    A2 = k2*(std::pow(2.0/beta, lambda) - std::pow(x0, lambda))/lambda;
  } else {
    k2=0.0;
    A2=0.0;
  }
  double k3 = std::pow(x_star, lambda-1.0);
  double A3 = 2.0*k3*std::exp(-x_star*beta/2.0)/beta;
  double A = A1 + A2 + A3;
  
  double U, V, X, h;
  do {
    U = R::runif(0.0, 1.0);
    V = R::runif(0.0, A);
    if(V <= A1){
      X = x0*V/A1;
      h = k1;
    } else if (V <= A1 + A2) {
      V -= A1;
      X = std::pow(std::pow(x0, lambda) + V*lambda/k2, 1.0/lambda);
      h = k2*std::pow(X, lambda-1);
    } else {
      V -= (A1+A2);
      X = -(2.0/beta) * std::log(std::exp(-x_star*beta/2.0) - V*beta/(2.0*k3));
      h = k3*std::exp(-X*beta/2);
    }
  } while(U*h > g(X, lambda, beta));
  
  return X;
}

double algo2(double const &lambda, double const &beta){
  double m = mode(lambda, beta);
  double temp = 1.0+lambda;
  double x_plus = (temp + std::sqrt(temp*temp + beta*beta)) / beta;
  double v_plus = std::sqrt(g(m, lambda, beta));
  double u_plus = x_plus*std::sqrt(g(x_plus, lambda, beta));
  
  double U, V, X;
  do {
    U = R::runif(0.0, u_plus);
    V = R::runif(0.0, v_plus);
    X = U/V;
  } while(V*V > g(X, lambda, beta));
  
  return X;
}


double algo3(double const &lambda, double const &beta){
  double xm, nc;
  double s, t; 
  double U, V, X; 
  
  double a, b, c;
  double p, q;
  double fi, fak;
  
  double y1, y2;
  
  double uplus, uminus; 
 
  t = 0.5 * (lambda-1.);
  s = 0.25 * beta;
  
  xm = mode(lambda, beta);
  
  nc = t*std::log(xm) - s*(xm + 1./xm);
  a = -(2.*(lambda+1.)/beta + xm);
  b = (2.*(lambda-1.)*xm/beta - 1.);
  c = xm;
  
  p = b - a*a/3.;
  q = (2.*a*a*a)/27. - (a*b)/3. + c;
  
  fi = std::acos(-q/(2.*sqrt(-(p*p*p)/27.)));
  fak = 2.*std::sqrt(-p/3.);
  y1 = fak * std::cos(fi/3.) - a/3.;
  y2 = fak * std::cos(fi/3. + 4./3.*M_PI) - a/3.;
  
  uplus  = (y1-xm) * std::exp(t*std::log(y1) - s*(y1 + 1./y1) - nc);
  uminus = (y2-xm) * std::exp(t*std::log(y2) - s*(y2 + 1./y2) - nc);
  
  do{
    U = R::runif(uminus, uplus);
    V = R::runif(0.0, 1.0);
    X = U/V + xm;
  } 
  while ((X <= 0.) || ((std::log(V)) > (t*std::log(X) - s*(X + 1./X) - nc)));
  
  
  return X;
}


double g(double const &x, double const &lambda, double const &beta){
  // quasi density
  return std::pow(x, lambda-1.0)*std::exp(-0.5*beta*(x+(1.0/x)));
}

double mode(double const &lambda, double const &beta){
  // mode of distribution
  if (lambda >= 1){
    return (sqrt((lambda-1.)*(lambda-1.) + beta*beta)+(lambda-1.))/beta;
  } else {
    return beta / (sqrt((1.-lambda)*(1.-lambda) + beta*beta)+(1.-lambda));
  }
}


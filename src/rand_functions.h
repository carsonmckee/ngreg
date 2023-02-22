#ifndef RAND_FUNCTIONS_H
#define RAND_FUNCTIONS_H

#include <RcppEigen.h>

Eigen::MatrixXd rmvnorm(const int &n, const Eigen::MatrixXd &mean, const Eigen::MatrixXd &cov);

double rgig1(double const &lamda, double const &chi, double const &psi);

double algo1(double const &lambda, double const &beta);

double algo2(double const &lambda, double const &beta);

double algo3(double const &lambda, double const &beta);

double g(double const &x, double const &lambda, double const &beta);

double mode(double const &lambda, double const &beta);

#endif
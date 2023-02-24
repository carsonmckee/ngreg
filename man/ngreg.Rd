\encoding{UTF-8}
\name{NG}
\alias{ngreg-package}
\alias{ngreg}
\alias{NG}
\title{
Bayesian Sparse Regression using Normal-Gamma Priors
}
\description{
  This package exposes a single method, \code{NG()}, which performs Bayesian linear regression using the Normal-Gamma priors described by Griffin & Brown (2010). This is particularly suited to problems where the number of covariates is greater than the number of observations. The MCMC sampler is implemented in C++ for performance, utilizing the template library 'Eigen' for fast linear algebra operations (interfaced to R using Rcpp/RcppEigen). In addition, OpenMP is enabled to allow Eigen to exploit multiple cores (on machines that have them) when performing dense matrix multiplication.
}

\usage{
NG(y, X, n_samples=10000, n_thin=1, init_alpha=0.0, init_beta=NULL, 
   init_psi=NULL, init_sigma=1.0, init_lambda=1.0, init_gamma=1.0,    
   init_prop_sd=1, tuning_time=2000, tuning_freq=200, verbose=TRUE)
}

\arguments{
  \item{y}{a numeric vector representing the response.}
  \item{X}{a numeric data.frame or matrix, representing the matrix of covariates. Should have the same number of rows as y. This should not contain an intercept column of ones as this will be added internally.}
  \item{n_samples}{the number of MCMC samples to be collected.}
  \item{n_thin}{the number of MCMC samples to discard before collecting the next sample. For example, if n_thin=4, then every 4th sample will be stored. By default, n_thin=1, corresonding to no thinning (all samples are collected).}
  \item{init_alpha}{the inital state of the markov chain for the intercept.}
  \item{init_beta}{the inital states of the markov chains for the regression coefficients. If supplied, it should be a numeric vector with length equal to the number of columns of X. If not supplied, then each chain is intialized at 0.}
  \item{init_psi}{the inital states of the markov chains for the latent psi variables. If supplied, it should be a numeric vector with length equal to the number of columns of X. If not supplied, then each chain is intialized to 1.}
  \item{init_sigma}{the inital state of the markov chain for the errors standard deviation. If not supplied, is set to 1.}
  \item{init_lambda}{the inital state of the markov chain for the lambda parameter. If not supplied, is set to 1.}
  \item{init_gamma}{the inital state of the markov chain for the gamma parameter. If not supplied, is set to 1.}
  \item{init_prop_sd}{the inital standard deviation for the Gaussian proposal distribution used in the Metropolis sampling step for lambda.}
  \item{tuning_time}{the number of initial iterations to be used for tuning the proposal standard deviation for the Metropolis step.}
  \item{tuning_freq}{the frequency with which to update the proposal standard deviation during the tuning_time.}
  \item{verbose}{if verbose=TRUE (default) then the percentage of MCMC samples collected will be displayed as sampling takes place, along with the time taken at the end. If FALSE, then nothing will be displayed.}
}

\value{
  \code{NG} returns an object of \link{class} "NG", which is essentially a holder containing the MCMC samples of the model parameters. 

An object of class "NG" contains the following components:
  \item{y}{the response used.}
  \item{X}{the model matrix used.}
  \item{call}{the matched call.}
  \item{n_samples}{the number of samples obtained.}
  \item{samples}{A dataframe holding the samples for each chain as columns.}
  \item{n_thin}{the level of thinning.}
  \item{tuning_time}{the number of iterations spent tuning the proposal variance for lambda.}
  \item{tuning_freq}{the tuning frequency used.}
  

}

\details{
  A summary of the model specification and sampling procedure is given below, however, full details can be found in Griffin & Brown (2010). The response vector, \eqn{\boldsymbol{Y}}, is modelled as: 
  
  \deqn{%
  \boldsymbol{Y} = \alpha + \boldsymbol{X}\boldsymbol{\beta} + \boldsymbol{\epsilon} \\
  }{}

  Where \eqn{\boldsymbol{X}} is an \eqn{n \times p} matrix of covariates, \eqn{\boldsymbol{\beta} = (\beta_1, ..., \beta_p)^T} is a vector of coefficients, \eqn{\boldsymbol{\epsilon} = (\epsilon_1, ..., \epsilon_n)^T}, with \eqn{\epsilon_i \overset{iid}{\sim} N(0, \sigma^2)} and \eqn{\alpha} is an intercept parameter. The hierarchical priors for these parameters are:
  
  \deqn{%
  \beta_j | \Psi_j \sim N(0, \Psi_j) \\
  \Psi_j | \lambda, \gamma \sim Ga(\lambda, \frac{1}{\gamma^2}) \\
  \lambda \sim Exp(1) \\
  \pi(\alpha) \propto 1 \\
  \pi(\sigma^{-2}) \propto 1
  }{}
  The prior for \eqn{\gamma} given \eqn{\lambda} is given by
  \deqn{%
  \nu_{\beta} = \lambda \gamma^2 \sim IG(2, M)
  }{}
  where \eqn{M = \frac{1}{p}\sum_{j=1}^{p} \hat{\beta_j}}. If X is non-singular (occurs when p > n-1), then \eqn{\hat{\beta}} is given by the least squares estimate, otherwise it is the Minimum Length Least Squares estimate.
  
  The MCMC sampler implements the Gibbs sampler as described in Griffin & Brown (2010), differing only in the adaptive Metropolis step for lambda. The adaptive step here uses the approach of Haario & Saksman (1999), where the proposal variance at the current iteration is based on the covariance of the previous \eqn{H} sampled values. This \eqn{H}, is denoted by \code{tuning_freq} here.
  
  The Gibbs sampling step for \eqn{\Psi_j} requires sampling from the Generalized Inverse Gaussian distribution.
  
  The function \code{summary} can be used to obtain a summary table containing posterior summaries, see \link{summary.NG} on how to use this.
}

\author{
Carson McKee \email{carson.mckee@kcl.ac.uk}
}

\keyword{Bayesian, Regression, Shrinkage, Sparse, Penalized}

\examples{
set.seed(1234)
p <- 100
n <- 75
X <- matrix(rnorm(p*n, 0, 1), nrow=n, ncol=p)
y <- 2 + 5*(X[,1] + X[,10] + X[,20] + X[,30] + X[,40]) + rnorm(n, 0, 1)

mod <- NG(y, X, verbose=FALSE)
plot(mod$samples$beta_1, type='l')

summary(mod)
}

\references{
Dirk Eddelbuettel, & Romain Francois (2011). Rcpp: Seamless R and C++ Integration. Journal of Statistical Software, 40(8), 1-18.

Douglas Bates, & Dirk Eddelbuettel (2013). Fast and Elegant Numerical Linear Algebra Using the RcppEigen Package. Journal of Statistical Software, 52(5), 1-24.

Gael Guennebaud, Benoit Jacob, & others. (2010). Eigen v3.

Haario, H., & Saksman, E. (1999). Adaptive Proposal Distribution for Random Walk Metropolis Algorithm. Computational Statistics, 14.

Hormann, J. (2014). Generating generalized inverse Gaussian random variates. Statistics and Computing, 24(4), 547-557.

Philip J. Brown, & Jim E. Griffin (2010). Inference with normal-gamma prior distributions in regression problems. Bayesian Analysis, 5(1), 171-188.
}

\seealso{
\code{\link{summary.NG}}, \code{\link{lm}}, \code{\link{glm}}
}

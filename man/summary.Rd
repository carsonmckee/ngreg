\encoding{UTF-8}
\name{summary.NG}
\alias{summary.NG}
\alias{summary}
\title{
Summarizing NG Model Objects
}
\description{
\code{summary} method for the class "NG". Calling this on an object of class 'NG' will return an object of class 'summary.NG'.
}

\usage{
\method{summary}{NG}(object, burnin=2000, estimate="mean", CI.level=0.05, nonzero_only=TRUE, ...)
}

\arguments{
  \item{object}{a fitted NG model object. See \link{NG}}
  \item{burnin}{the number of initial MCMC iterations to discard as burn-in when calculating posterior estimates.}
  \item{estimate}{the posterior estimate to use. Valid values are "mean" and "median".}
  \item{CI.level}{the level of the credible intervals reported by the summary object.}
  \item{nonzero_only}{if TRUE then only the regression coefficients whose credible interval does not contain zero are displayed.}
  \item{...}{further arguments passed to or from other methods.}
}

\value{
  An object of class 'summary.NG' containing the following:
  \item{call}{the matched call.}
  \item{nonzero_only}{a boolean value indicating whether to only store regression coefficients whose credible intervals do not contains zero.}
  \item{coefficients}{a table containing the posterior estimate (mean or median) and credible intervals (at the CI.level) for the regression coefficients, lambda, gamma and sigma (error standard deviation) parameters.}
}

\details{
  This method provides posterior summaries for the regression coefficients and other model parameters. It does not include summaries for the latent \eqn{\Psi_j} parameters.
}

\author{
Carson McKee \email{carson.mckee@kcl.ac.uk}
}


\keyword{Bayesian, Regression, Shrinkage, Sparse, Penalized}


\examples{
p <- 100
n <- 75
X <- matrix(rnorm(p*n, 0, 1), nrow=n, ncol=p)
y <- 2 + 5*(X[,1] + X[,10] + X[,20] + X[,30] + X[,40]) + rnorm(n, 0, 1)

mod <- NG(y, X, verbose=FALSE)
summary(mod, estimate='median', CI.level=0.01, nonzero_only=TRUE)
}


\seealso{
  \code{\link{NG}}
}


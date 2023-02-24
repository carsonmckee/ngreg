
summary.NG <- function(object, burnin=2000, estimate = "mean", CI.level =
                         0.05, nonzero_only = TRUE, ...){
  
  if(estimate == "mean"){
    est <- mean
  } else if(estimate == "median"){
    est <- median
  } else {
    stop('Only mean and median are supported for posterior estimates.')
  }
  
  if(CI.level <= 0 | CI.level >= 1){
    stop('CI.level must be in range (0, 1)')
  }
  
  p <- ncol(object$X)
  chains <- object$samples[burnin:object$n_samples,-((p+2):(2*p+1))]
  coefficients <- cbind(apply(chains, 2, est),
                        t(apply(chains, 
                                2, 
                                quantile, 
                                probs=c(CI.level/2, 1-CI.level/2)))
                        )
  colnames(coefficients)[1] <- estimate
  info <- ''
  if(nonzero_only){
    coefficients <- coefficients[!((coefficients[,2] < 0) & (coefficients[,3] > 0)),]
    info <- 'Showing non zero coefficients only.'
  }
  
  obj <- list(
    call = object$call,
    nonzero_only = nonzero_only,
    coefficients = coefficients 
  )
  attr(obj, 'class') <- 'summary.NG'
  
  obj
}

print.summary.NG <- function(x, ...){
  cat("\nCall:\n")
  print(x$call)
  if(x$nonzero_only){
    cat("\nCoefficients (non-zero only):\n")
  } else {
    cat("\nCoefficients:\n")
  }
  print(head(x$coefficients, nrow(x$coefficients)-3))
  cat("\nErrors Standard Deviation:\n")
  print(x$coefficients['sigma',])
  cat("\nShrinkage Parameters:\n")
  print(tail(x$coefficients, 2))
  
}

print.NG <- function(x, ...){
  cat("\nCall:\n")
  print(x$call)
  cat("\nModel Info:\n")
  cat("n         =", nrow(x$X), "\n")
  cat("p         =", ncol(x$X), "\n")
  cat("\nMCMC Info:\n")
  cat("n_samples   =", x$n_samples, "\n")
  cat("n_thin      =", x$n_thin, "\n")
  cat("tuning_time =", x$tuning_time, "\n")
  cat("tuning_freq =", x$tuning_freq, "\n")
  cat("*** Call summary(NG_OBJECT) to obtain posterior summaries ***")
}

NG <- function(
    y, 
    X, 
    n_samples    = 10000,
    n_thin       = 1,
    init_alpha   = 0.0, 
    init_beta    = NULL,
    init_psi     = NULL, 
    init_sigma   = 1.0, 
    init_lambda  = 1.0,
    init_gamma   = 1.0, 
    init_prop_sd = 1,
    tuning_time  = 2000, 
    tuning_freq  = 200, 
    verbose      = TRUE){
  
  if(verbose){
    tictoc::tic()
  }
  
  call_ <- match.call()
  
  p <- ncol(X)
  n <- nrow(X)
  
  if(is.data.frame(X)){
    col_names <- colnames(X)
    psi_names <- paste(rep('psi_', p), col_names, sep='')
    x <- as.matrix(X)
  } else if(is.matrix(X)) {
    col_names <- paste(rep('beta_', p), 1:p, sep='')
    psi_names <- paste(rep('psi_', p), 1:p, sep='')
    x <- X
  } else {
    stop('x must be a dataframe or a matrix')
  }
  
  if(length(y) != n){
    stop('Response, y, has differing number of rows compared to design matrix.')
  }
  
  if(n_thin < 1){
    stop('n_thin must be a positive integer.')
  }
  
  if(p < n){
    beta_hat <- lm(y ~ x)$coefficients[2:(p+1)]
  } else {
    beta_hat <- qr.coef(qr(x), y)
    beta_hat[is.na(beta_hat)] <- 0
  }
  M <- mean(beta_hat^2)
  x_with_intercept <- cbind(rep(1, n), x) 
  X_t_X <- t(x_with_intercept) %*% x_with_intercept # only need to do this once
  X_t_y <- t(x_with_intercept) %*% y # only need to do this once
  
  # if no initial positions supplied for beta/psi then set all to 0/1.
  if(is.null(init_beta)){
    init_beta <- matrix(rep(0, p), nrow=p, ncol=1)
  }
  if(is.null(init_psi)){
    init_psi <- matrix(rep(1, p), nrow=1, ncol=p)
  }
  
  chains <- .Call('_ngreg_get_chains', PACKAGE = 'ngreg', 
                  y = matrix(y, n, 1),
                  x = x,
                  X_t_X = X_t_X,
                  X_t_X = X_t_y,
                  n_samp = n_samples,
                  n_thin = n_thin,
                  n = n,
                  p = p,
                  M = M,
                  init_alpha = init_alpha,
                  init_beta = init_beta,
                  init_psi = init_psi,
                  init_sigma = init_sigma,
                  init_lambda = init_lambda,
                  init_gamma = init_gamma, 
                  init_prop_sd=init_prop_sd,
                  tuning_time = tuning_time,
                  tuning_freq = tuning_freq,
                  verbose = as.integer(verbose))
  
  chains <- data.frame(chains)
  names(chains) <- c('alpha', 
                     col_names, 
                     psi_names, 
                     'sigma',
                     'lambda',
                     'gamma')
  
  # create NG object
  obj <- list(
    call=call_,
    y = y,
    X = X,
    n_samples=n_samples,
    n_thin=n_thin,
    tuning_freq=tuning_freq,
    tuning_time=tuning_time,
    samples=chains
  )
  attr(obj, 'class') <- 'NG'
  
  gc()
  
  if(verbose){
    tictoc::toc()
  }
  
  obj
}

# set.seed(1234)
# p <- 100
# n <- 10000
# X <- matrix(rnorm(p*n, 0, 1), nrow=n, ncol=p)
# y <- 2 + 5*(X[,1] + X[,10] + X[,20] + X[,30] + X[,40]) + rnorm(n)
# 
# mod <- NG(y, X, n_samples=50000)
# 
# summary(mod)

# plot(mod$samples$alpha, type='l')
# plot(mod$samples$beta_1, type='l')
# plot(mod$samples$gamma, type='l')
# plot(mod$samples$lambda, type='l')
# # plot(mod$samples$sigma, type='l')


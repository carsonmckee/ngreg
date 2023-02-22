
summary.NG <- function(obj, burn_in=2000, ...){
  # custom summary for NG object
  out <- list(
    call=obj$call
  )
}


print.NG <- function(object, ...){
  cat("\nCall:\n")
  print(object$call)
}

create_value <- function(chains, 
                         n_samples, 
                         n_thin, 
                         tuning_freq, 
                         tuning_time,
                         call_=NULL){
  obj <- list(
    call=call_,
    n_samples=n_samples,
    n_thin=n_thin,
    tuning_freq=tuning_freq,
    tuning_time=tuning_time,
    samples=chains,
    CI = t(apply(chains , 2 , quantile , probs=c(0.05, 0.95)))
  )
  attr(obj, 'class') <- 'NG'
  
  obj
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
    init_prop_sd = 0.1,
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
  
  obj <- create_value(chains=chains, 
                      n_samples=n_samples, 
                      n_thin=n_thin, 
                      tuning_freq=tuning_freq, 
                      tuning_time=tuning_time,
                      call=call_)
  gc()
  
  if(verbose){
    tictoc::toc()
  }
  return(obj)
}

# set.seed(1234)
# p <- 80
# n <- 200
# X <- matrix(rnorm(p*n, 0, 1), nrow=n, ncol=p)
# y <- -2 + 1*(X[,1] + X[,10] + X[,30] + X[,40]) + rnorm(n, 0, 1)
# 
# mod <- NG(y, X)
# 
# plot(mod$samples$alpha, type='l')
# plot(mod$samples$beta_1, type='l')
# plot(mod$samples$beta_2, type='l')
# plot(mod$samples$beta_10, type='l')
# plot(mod$samples$sigma, type='l')
# plot(mod$samples$gamma, type='l')
# plot(mod$samples$lambda, type='l')
# 
# summary(mod)

# mod$CI[!((mod$CI[,1] < 0) & (mod$CI[,2] > 0)),]
# acf(get_chains(y, X, 10000, n_thin=1, verbose=TRUE)$beta_1)
# acf(get_chains(y, X, 10000, n_thin=10, verbose=TRUE)$beta_1)


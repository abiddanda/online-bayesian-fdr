
#' Gibbs sampler to estimate mixture model 
#' @param X vector of test.statistics at time t
#' @param n_iter number of iterations of the gibbs sampler
#' @param n_burnin number of burn in iterations
#' @param sampling_interval interval of sampling times
#' @return 
gibbs_sampler <- function(X, n_iter, n_burnin=50, sampling_interval=10){
  n <- length(X)
  
  # intialize params
  a <- 1 # gamma prior hyper param
  b <- 1 # gamma prior hyper param
  m <- 0 # normal prior hyper param
  alpha <- 1
  beta <- 1
  alpha_norm <- 2 # normal hyper param 
  pi0 <- rbeta(1, alpha, beta)
  phi1 <- rgamma(1, a/2, b/2)
  mu1 <- rnorm(1, m, 1 / (alpha_norm * phi1))
  
  pi0_samples <- c()
  phi1_samples <- c()
  mu1_samples <- c()
  
  trace <- rep(NA, n_iter)
  
  # samples
  for(i in 1:n_iter){
    # update z 
    p_z_num <- pi0 * exp(-.5*(X^2))
    p_z_denom <- p_z_num + ((1 - pi0) * phi1 * exp(-phi1/2 * (X - mu1)^2))
    p_z0 <- p_z_num / p_z_denom
    z <- rbinom(n, 1, 1-p_z0)
    
    # update pi0
    pi0 <- rbeta(1, alpha + sum(z == 0), beta + sum(z == 1))
    
    # print(sum((X[z==1] - mu1)^2))
    # update phi1
    phi1 <- rgamma(1, (a + sum(z == 1)) / 2, (b + sum((X[z==1] - mu1)^2)) / 2)
    
    # update mu1
    m <- ((alpha_norm * m) + sum(X[z==1])) / (alpha_norm + sum(z == 1))
    mu1 <- rnorm(1, m, 1 / ((alpha_norm + sum(z == 1)) * phi1))
    
    # Calculating the trace
    trace[i] <- sum(log(p_z0[which(z==0)])) + sum(log((1-p_z0)[which(z==1)]))

    if (i > n_burnin){
      pi0_samples <- c(pi0_samples, pi0)
      phi1_samples <- c(phi1_samples, phi1)
      mu1_samples <- c(mu1_samples, mu1)
    }
  }
  return(list(pi0_samples = pi0_samples, phi1_samples = phi1_samples, 
              mu1_samples = mu1_samples, trace=trace))
}



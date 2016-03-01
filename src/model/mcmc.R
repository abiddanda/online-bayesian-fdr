library(ggplot2)
source("../simulation/simulate_mixture.R")


#' Gibbs sampler to estimate mixture model 
#' @param X vector of test.statistics at time t
#' @param n_iter number of iterations of the gibbs sampler
#' @param n_burnin number of burn in iterations
#' @param sampling_interval interval of sampling times
#' @return 
gibbs_sampler <- function(X, n_iter, n_burnin, sampling_interval){
  n <- length(X)
  
  # intialize params
  a <- 1 # gamma prior hyper param
  b <- 1 # gamma prior hyper param
  m <- 1 # normal prior hyper param
  alpha <- 1
  beta <- 1
  alpha_norm <- 2 # normal hyper param 
  pi0 <- rbeta(1, alpha, beta)
  phi1 <- rgamma(1, a/2, b/2)
  mu1 <- rnorm(1, m, 1 / (alpha_norm * phi1))
  
  pi0_samples <- c()
  phi1_samples <- c()
  mu1_samples <- c()
  
  # sammples
  for(i in 1:n_iter){
    # update z 
    p_z_num <- pi0 * exp(-.5*(X^2))
    p_z_denom <- p_z_num + ((1 - pi0) * phi1 * exp(-phi1/2 * (X - mu1)^2))
    p_z <- p_z_num / p_z_denom
    z <- rbinom(n, 1, p_z)
    
    # update pi0
    pi0 <- rbeta(1, alpha + sum(z == 0), beta + sum(z == 1))
    
    # update phi1
    phi1 <- rgamma(1, (a + sum(z == 1)) / 2, (b + sum((X[z==1] - mu1)^2)) / 2)
    
    # update mu1
    m <- ((alpha_norm * m) + sum(X[z==1])) / (alpha_norm + sum(z == 1))
    mu1 <- rnorm(1, m, 1 / ((alpha_norm + sum(z == 1)) * phi1))
    
    if (i > n_burnin){
      print(i)
      pi0_samples <- c(pi0_samples, pi0)
      phi1_samples <- c(phi1_samples, phi1)
      mu1_samples <- c(mu1_samples, mu1)
    }
  }
  return(list(pi0_samples = pi0_samples, phi1_samples = phi1_samples, 
              mu1_samples = mu1_samples))
}

### RUN STUFF
#simulation_list <- sim.mixture.2comp(n = 1000, pi0 = .9, mu1 = 2, sigma1 = 1)
#X <- simulation_list$Z
#gibbs_list <- gibbs_sampler(X, n_iter = 5000, n_burnin = 100, sampling_interval = 10)
#plot(gibbs_list$mu1_samples, type="l")
#hist(gibbs_list$mu1_samples, breaks=10)
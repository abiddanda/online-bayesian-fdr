# ----------------------------------
# Different FDR Methods in R
# ----------------------------------

# 1. Bayesian FDR method
# ----------------------------------
BayesFDR <- function(z, mu, alpha=0.2, pi0=0.8, sigma2=4){
  numerator <- pi0*(1 - pnorm(z))
  denom <- pi0*(1 - pnorm(z)) 
  denom <- denom + (1-pi0)*(1 - pnorm((z-mu)/sqrt(sigma2)))
  return(numerator/denom - alpha)
}

# 2. BH-FDR
# -----------------------------------
BHFDR <- function(z,mu, alpha=0.2, pi0=0.8, sigma2=4){
  numerator <- (1 - pnorm(z))
  denom <- pi0*(1 - pnorm(z)) 
  denom <- denom + (1-pi0)*(1 - pnorm((z-mu)/sqrt(sigma2)))
  return(numerator/denom - alpha)
}

# 3. Storey FDR
# -----------------------------------
StoreyFDR <- function(z, mu, alpha=0.2, pi0=0.8, sigma2=4){
  pi_hat <- pi0 + 2*(1-pi0)*pnorm(-mu/sqrt(sigma2))
  numerator <- (1 - pnorm(z)) * pi_hat
  denom <- pi0*(1 - pnorm(z)) 
  denom <- denom + (1-pi0)*(1 - pnorm((z-mu)/sqrt(sigma2)))
  return(numerator/denom - alpha)
}

# Plotting functions
# -----------------------------------

get_roots <- function(f, alpha, mus, pi0s, sigma2s){
  roots <- rep(NA, length(mus))
  for (i in 1:length(mus)){
    roots[i] <- uniroot(f, interval=c(0,8), extendInt = "yes", alpha=alpha, mu=mus[i], pi0=pi0s[i], sigma2=sigma2s[i])$root 
  }
  return(roots)
}

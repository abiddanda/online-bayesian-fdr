
# ------------------------------------
# Calculating the FPR and Power across
# ------------------------------------

# Sourcing Relevant Methods
source("fdr_methods.R")
source("lbond_lbord_ainvesting.R")


# Sourcing Data
X <- readRDS("~/Repos/online-bayesian-fdr/data/sim1.rds")

# Converting Z-scores to Pvals
convert_Z_to_Pval <- function(mixture.sim, alpha=0.1){
  Zs <- mixture.sim$Z
  return(pnorm(-abs(Zs)))
}

compute_empirical_FDP_LBOND <- function(mixture.sim, alpha=0.1){
  pvals <- convert_Z_to_Pval(mixture.sim)
  return(alpha_inv(pvals, alpha))
}

compute_empirical_FDP_LBORD <- function(mixture.sim, alpha=0.1){
  pvals <- convert_Z_to_Pval(mixture.sim)  
  return(LBOND(pvals, alpha))
}

compute_empirical_FDP_Ainvesting <- function(mixture.sim){
  
}


bayesian_Power <- function(pi0=0.9, mu=2, sigma2=1, alpha=0.1){
  z.hat <- get_roots(BayesFDR, alpha, mus=mu, pi0s=pi0, sigma2s=sigma2)
  prob_signals <- pi0 *(1-pnorm(z.hat)) + (1-pi0)*(1-pnorm((z.hat-mu)/sqrt(sigma2)))
  fdr <- BayesFDR(z.hat, mu, alpha=alpha, pi0 = pi0, sigma2 = 1) + alpha
  power <- (prob_signals * (1-fdr)) / (1 - pi0)
  print(sprintf("FDR : %f , Power : %f, Z.hat : %f", fdr, power, z.hat))
  # return(power)
}

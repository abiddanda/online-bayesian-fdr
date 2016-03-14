
# ------------------------------------
# Calculating the FPR and Power across
# ------------------------------------

# Sourcing Relevant Methods
source("~/Repos/online-bayesian-fdr/src/analysis/fdr_methods.R")
source("~/Repos/online-bayesian-fdr/src/analysis/lbond_lbord_ainvesting.R")
source("~/Repos/online-bayesian-fdr/src/model/mcmc.R")

# Reading in Data
sim10000 <- readRDS("~/Repos/online-bayesian-fdr/data/sim1.rds")
sim500 <- readRDS("~/Repos/online-bayesian-fdr/data/sim2_500_08_3_1.rds")

# Method for Sequential Gibbs Updates
sequential_gibbs <- function(data, timestep, n_iter=1000, n_burnin=50, sampling_interval=10){
  n <- length(data)
  if (n %% timestep != 0){stop("Not an even divisor!")}
  
  ncolumns <- (n_iter - n_burnin) / sampling_interval
  print(ncolumns)
  mu.mat <- matrix(NA, nrow=(n/timestep), ncol = ncolumns)
  pi0.mat <- matrix(NA, nrow=(n/timestep), ncol= ncolumns)
  sigma2.mat <- matrix(NA, nrow=(n/timestep), ncol=ncolumns)
  mat.ind <- 1
  for (t in seq(1, n, timestep)){
    print(t)
    data.t <- data[1:t]
    cur.param.samples <- gibbs_sampler(data.t, n_iter, sampling_interval = 10)
    mu.mat[mat.ind, ] <- cur.param.samples$mu1_samples
    pi0.mat[mat.ind, ] <- cur.param.samples$pi0_samples
    sigma2.mat[mat.ind, ] <- 1/ cur.param.samples$phi1_samples
    mat.ind <- mat.ind + 1  
  }
  return(list(mu.samples=mu.mat, pi0.samples=pi0.mat, sigma2.samples=sigma2.mat))
}

# Computing credible intervals of the different parameters
cred_intervals <- function(sum_stat_mat){
  n <- nrow(sum_stat_mat)
  cred.int <- matrix(NA, nrow=n, ncol=3)
  for (i in 1:n){
    cred.int[i, ] <- quantile(sum_stat_mat[i,], c(0.05,0.5,0.95))
    if (cred.int[i,1] < 0){cred.int[i,1] <- 0.0}
  }
  return(cred.int)
}


# FDP and Power for Current Methods
compute_empirical_FDP_LBOND <- function(mixture.sim, alpha=0.1){
  lbond_results <- LBOND(mixture.sim$Pval, alpha)
  lbond_power <- sum(((lbond_results$discoveries == 1) & (mixture.sim$true_signals == 1))) / sum((mixture.sim$true_signals == 1))
  lbond_fdp <- sum((lbond_results$discoveries == 1) & (mixture.sim$true_signals == 0)) / sum((mixture.sim$true_signals == 0))
  lbond_results$power <- lbond_power
  lbond_results$fdp <- lbond_fdp
  return(lbond_results)
}

compute_empirical_FDP_LBORD <- function(mixture.sim, alpha=0.1){
  lbord_results <- LBORD(mixture.sim$Pval, alpha)
  lbord_power <- sum(((lbord_results$discoveries == 1) & (mixture.sim$true_signals == 1))) / sum((mixture.sim$true_signals == 1))
  lbord_fdp <- sum((lbord_results$discoveries == 1) & (mixture.sim$true_signals == 0)) / sum((mixture.sim$true_signals == 0))
  lbord_results$power <- lbord_power
  lbord_results$fdp <- lbord_fdp
  return(lbord_results)
}

compute_empirical_FDP_Ainvesting <- function(mixture.sim, alpha=0.1){
  alpha_inv_results <- alpha_inv(mixture.sim$Pval, alpha)
  alpha_inv_power <- sum(((alpha_inv_results$discoveries == 1) & (mixture.sim$true_signals == 1))) / sum((mixture.sim$true_signals == 1))
  alpha_inv_fdp <- sum((alpha_inv_results$discoveries == 1) & (mixture.sim$true_signals == 0)) / sum((mixture.sim$true_signals == 0))
  alpha_inv_results$power <- alpha_inv_power
  alpha_inv_results$fdp <- alpha_inv_fdp
  return(alpha_inv_results)
}

# Getting the FDP and Power for our method
bayesian_criteria <- function(pi0s, mus, sigma2s, alpha=0.1){
  n <- length(pi0s)
  if (length(pi0s) != length(mus)){stop("Not good!")}
  z.hat <- rep(NA, n)
  for (i in 1:n){
    z.hat[i] <- get_roots(BayesFDR, alpha, mus=mu[i], pi0s=pi0s[i], sigma2s=sigma2s[i])
  }
  return(z.hat)
}

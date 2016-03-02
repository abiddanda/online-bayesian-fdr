#-----------------------------------
# Description : 
#   Simulates data from a two 
#   component mixture.
#
# Args:
#		pi0 - proportion from component1
#		mu0 - mean of nulls
#   sd1 - var of nulls
#   mu1 - mean of signals
#   sd1 - var of signals
#-----------------------------------
sim.mixture.2comp <- function(n, pi0, mu0=0, sigma0=1, mu1, sigma1){
	signals <- rbinom(n, 1, (1-pi0))
	z <- vapply(signals,
	            FUN = function(x){ifelse(x,rnorm(1,mu1,sigma1), rnorm(1,mu0,sigma0))}, 
	            FUN.VALUE = 1)
	pvals <- 2*pnorm(-abs((z-mu0)/sqrt(sigma0)))
	mixture.output <- list(true_signals = signals, Z=z, Pval=pvals)
	return(mixture.output)
}


#-----------------------------------
# Description : 
#   Simulates data from a two 
#   component mixture with correlation
#   between the transitions
#
# Args:
#		pi0 - proportion from component1
#		mu0 - mean of nulls
#   sd1 - var of nulls
#   mu1 - mean of signals
#   sd1 - var of signals
#-----------------------------------
sim.mixture.2comp.corr <- function(n, trans.mat, init.val=1, mu0=0, sigma0=1, mu1, sigma1){
  #1. Simulate signals from Markov Chain
  signals <- rep(NA, n)
  signals[1] <- init.val+1
  for (t in 2:n){
    p <- trans.mat[signals[t-1], ]
    signals[t] <- which(rmultinom(1, 1, p) == 1) 
  }
  # because R is 1-indexed
  signals <- signals-1
  z <- vapply(signals,
              FUN = function(x){ifelse(x,rnorm(1,mu1,sigma1), rnorm(1,mu0,sigma0))},
              FUN.VALUE = 1)
  pvals <- 2*pnorm(-abs((z-mu0)/sqrt(sigma0)))
  mixture.output <- list(true_signals = signals, Z=z, Pval=pvals)
  return(mixture.output)
}


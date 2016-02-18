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
	mixture.df <- data.frame(true_signals = signals, Z=z, P=pvals)
	return(mixture.df)
}


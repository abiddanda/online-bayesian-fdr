#-----------------------------------
# Description : 
#   Simulates data from a two 
#   component mixture.
#
# Args:
#		pi0 - proportion from component1
#		mu0 - mean of nulls
#   sd1 - sd of nulls
#   mu1 - mean of signals
#   sd1 - sd of signals
#-----------------------------------
sim.mixture.2comp <- function(n, pi0, mu0, sd0, mu1, sd1){
	signals <- rbinom(n, 1, (1-pi0))
	x <- vapply(signals,
	            FUN = function(x){ifelse(x,rnorm(1,mu1,sd1), rnorm(1,mu0,sd0))}, 
	            FUN.VALUE = 1)
	# TODO : simulate p-values from the mixture?
	pvals <- 2*pnorm(-abs((x-mu0)/sd0))
	mixture.df <- data.frame(true_signals = signals, X=x, P=pvals)
	return(mixture.df)
}


# TODO :
#  - grouped signals
















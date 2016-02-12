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
	x <- vapply(rbinom(n, 1, (1-pi0)),
	            FUN = function(x){ifelse(x,rnorm(1,mu0,sd0), rnorm(1,mu1,sd1))}, 
	            FUN.VALUE = 1)
	return(x)
}

# TODO :
#  - 3 component mixtures
#  - grouped signals
















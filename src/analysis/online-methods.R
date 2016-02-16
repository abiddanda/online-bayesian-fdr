
#--------------------------- 
# Estimate pi0 (ala Storey)
# Args:
#   p = current vector of 
#     p-values
#   lambda = lambda value 
#        from Storey 2002
#---------------------------
estimate.nulls <- function(p, lambda = 0.5){
  pi0 <- sum(p > lambda) / ((1-lambda)*length(p))
  return(pi0)
}


# --------------------------
# Sequentially Estimate pi0
# Args:
# 	p = full vector of p-values
# 	lambda = lambda value from
# 			Storey 2002
# 	tau = min. starting point
# --------------------------
estimate.nulls.sequential <- function(p, lambda=0.5, tau=100){
  if (length(p) < tau){stop()}
  n <- length(p)
  pi0.t <- vapply(seq(tau,n), FUN = function(x){estimate.nulls(p[1:x], lambda=lambda)}, FUN.VALUE = 1)
  return(pi0.t)
}


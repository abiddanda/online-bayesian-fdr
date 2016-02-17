
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
estimate.nulls.seq <- function(p, lambda=0.5, tau=100){
  if (length(p) < tau){stop("Not enough p-values!")}
  n <- length(p)
  pi0.t <- vapply(seq(tau,n), FUN = function(x){estimate.nulls(p[1:x], lambda=lambda)}, FUN.VALUE = 1)
  return(pi0.t)
}


# --------------------------
# Estimate pFDR (from Storey)
# Args:
#   p = vector of sequential 
#     p-values
#   pi0 = sequential prop. 
#     nullss
# --------------------------
estimate.pfdr <- function(p, pi0, gamma){
  m <- length(p)
  R <- sum(p < gamma)
  pr.gamma <- ifelse(R == 0, 1, R) / m
  return((pi0*gamma)/ (pr.gamma * (1 - (1-gamma)^m)))
}

estimate.pfdr.seq <- function(p, pi0.t, gamma, tau=100){
  n <- length(p)
  pfdr.t <- vapply(seq(tau,n), FUN = function(x){estimate.pfdr(p[1:x], pi0.t[x-tau+1], gamma)}, FUN.VALUE = 1)
  return(pfdr.t)
}

# --------------------------
# Estimate FDR (from Storey)
# Args:
#   p = vector of p-values
#   pi0 = current estimate 
#     of pi0
# --------------------------
estimate.fdr <- function(p, pi0, gamma){
  m <- length(p)
  R <- sum(p < gamma)
  pr.gamma <- ifelse(R==0, 1, R) / m
  return((pi0*gamma)/(pr.gamma))
}

estimate.fdr.seq <- function(p, pi0.t, gamma, tau=100){
  n <- length(p)
  pfdr.t <- vapply(seq(tau,n), FUN = function(x){estimate.fdr(p[1:x], pi0.t[x-tau+1], gamma)}, FUN.VALUE = 1)
  return(pfdr.t)
}

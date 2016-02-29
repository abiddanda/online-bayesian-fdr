alpha_inv = function(P,alpha){
	# alpha-investing by Foster & Stine (2007)
	# W(0) = alpha
	# at time i reject if P[i]<=alpha[i]
	# choose alpha[i]<=W(i-1)/(1+W(i-1))
	# if reject, W(i) = W(i-1) + alpha
	# if not, W(i) = W(i-1) - alpha[i]/(1-alpha[i])
	
	# scheme: choose alpha[i] = 0.5*W(i-1)/(1+W(i-1))
	
	n = length(P)
	W = rep(0,n+1)
	alphas = rep(0,n)
	W[1] = alpha
	discoveries = rep(0,n)
	
	for(i in 1:n){
		alphas[i] = 0.5 * W[i] / (1+W[i])
		if(P[i]<=alphas[i]){
			discoveries[i] = 1
			W[i+1] = W[i]+alpha
		}else{
			W[i+1] = W[i] - alphas[i]/(1-alphas[i])
		}		
	}
	output = list()
	output$discoveries = discoveries
	output$alphas = alphas
	output
}

LBOND = function(P,alpha){
	# "level based on number of discoveries", Javanmard & Montanari 2015
	# at time i reject if P[i]<=alpha[i]
	# set alpha[i] = beta[i]*max{1,D(i-1)}
	# where D(i-1) is the number of discoveries up to time i-1
	# and beta is a sequence with sum_{i=1}^{infty} beta[i] = alpha
	
	# scheme: beta[i] proportional to i^{-1.5}
	
	n = length(P)
	alphas = rep(0,n)
	beta = (1:n)^(-1.5); beta = beta * alpha / sum(beta)
	discoveries = rep(0,n)
	ndisc = 0
	
	for(i in 1:n){
		alphas[i] = beta[i] * max(1,ndisc)
		if(P[i]<=alphas[i]){
			discoveries[i] = 1
			ndisc = ndisc + 1
		}
	}
	
	output = list()
	output$discoveries = discoveries
	output$alphas = alphas
	output	
}

LBORD = function(P,alpha){
	# "level based on recent discoveries", Javanmard & Montanari 2015
	# at time i reject if P[i]<=alpha[i]
	# set alpha[i] = beta[i-tau(i)]
	# where tau(i) is the time of the most recent discovery
	# and beta is a sequence with sum_{i=1}^{infty} beta[i] = alpha
	
	# scheme: beta[i] proportional to i^{-1.5}
	
	n = length(P)
	alphas = rep(0,n)
	beta = (1:n)^(-1.5); beta = beta * alpha / sum(beta)
	discoveries = rep(0,n)
	tau = 0
	
	for(i in 1:n){
		alphas[i] = beta[i - tau]
		if(P[i]<=alphas[i]){
			discoveries[i] = 1
			tau = i
		}
	}
	
	output = list()
	output$discoveries = discoveries
	output$alphas = alphas
	output	
}



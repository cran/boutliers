STR <- function(y, v, method="REML", data=NULL, B=2000, alpha=0.95, seed=123456){

	set.seed(seed)

	if(is.null(data)==FALSE){

		data <- data.frame(data)

		y <- data[, deparse(substitute(y))]
		v <- data[, deparse(substitute(v))]

	}

	reml1 <- rma(yi=y,vi=v,method=method)

	mu0 <- reml1$beta
	v0 <- v + reml1$tau2

	n <- length(y)
	psi <- numeric(n)

	for(i in 1:n){
	
		y_i <- y[setdiff(1:n,i)]
		v_i <- v[setdiff(1:n,i)]
	
		reml_i <- rma(yi=y_i,vi=v_i,method=method)
		
		W0_i <- (v + reml_i$tau2)^-1
		v_psi_i <- (W0_i[i])^-1 + sum(W0_i[-i])^-1
	
		psi[i] <- (y[i] - reml_i$beta) / sqrt(v_psi_i)
		
	}
                              
	psi.b <- matrix(numeric(n*B),B)

	for(i in 1:n){
	
		y_i <- y[setdiff(1:n,i)]
		v_i <- v[setdiff(1:n,i)]
	
		reml_i <- rma(yi=y_i,vi=v_i,method=method)

		for(b in 1:B){
	
			y.b <- rnorm(n, mean=reml_i$beta, sd=sqrt(v+reml_i$tau2))
		
			y.b_i <- y.b[setdiff(1:n,i)]
	
			reml.b_i <- rma(yi=y.b_i,vi=v_i,method=method)
		
			W0_i <- (v + reml_i$tau2)^-1
			v_psi_i <- (W0_i[i])^-1 + sum(W0_i[-i])^-1
	
			psi.b[b,i] <- (y.b[i] - reml_i$beta) / sqrt(v_psi_i)

		}
	
		print1 <- paste0("The bootstrap for ",i,"th study is completed.")
		if(b%%100==0) message(print1)

	}
		
	Q1 <- Q2 <- numeric(n)
	
	alpha1 <- (1 - alpha)/2
	alpha2 <- 1 - alpha1
		
	for(i in 1:n){
	
		X.b <- psi.b[,i]
		#P[i] <- QT(X.b, abs(psi[i]))
		Q1[i] <- as.numeric(quantile(X.b,alpha1))
		Q2[i] <- as.numeric(quantile(X.b,alpha2))
	
	}

	id <- 1:n
	R <- data.frame(id,psi,Q1,Q2)
	R <- R[rev(order(abs(psi))),]

	return(R)

}


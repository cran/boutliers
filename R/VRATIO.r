VRATIO <- function(y, v, B=2000, alpha=0.05){

	reml1 <- rma(yi=y,vi=v)
	
	mu0 <- as.numeric(reml1$beta)
	v0 <- v + reml1$tau2

	V0 <- reml1$tau2
	V1 <- reml1$se^2

	#reml1 <- REML(y,v)
	
	#mu0 <- reml1$mu
	#v0 <- v + reml1$V0

	#V0 <- reml1$V0
	#V1 <- reml1$V1

	n <- length(y)
	VR <- TR <- numeric(n)

	for(i in 1:n){
	
		y_i <- y[setdiff(1:n,i)]
		v_i <- v[setdiff(1:n,i)]
	
		#reml_i <- REML(y_i,v_i)
		
		#VR[i] <- reml_i$V1 / V1
		#TR[i] <- reml_i$V0 / V0
		
		reml_i <- rma(yi=y_i,vi=v_i)
		
		VR[i] <- reml_i$se^2 / V1
		TR[i] <- reml_i$tau2 / V0		
		
	}
                              
	VR.b <- TR.b <- matrix(numeric(n*B),B)

	for(b in 1:B){
	
		y.b <- rnorm(n, mean=mu0, sd=sqrt(v0))
		reml.b <- REML(y.b,v)
		
		V0.b <- reml.b$V0
		V1.b <- reml.b$V1
		
		for(i in 1:n){
	
			y_i <- y.b[setdiff(1:n,i)]
			v_i <- v[setdiff(1:n,i)]
	
			reml_i <- REML(y_i,v_i)
		
			VR.b[b,i] <- reml_i$V1 / V1.b
			TR.b[b,i] <- reml_i$V0 / V0.b

		}
	
		print1 <- paste0("The ",b,"th bootstrap is completed.")
		if(b%%100==0) print(print1)

	}
		
	Q1 <- Q2 <- numeric(n)	
		
	for(i in 1:n){
	
		X1 <- VR.b[,i]
		X2 <- TR.b[,i]
		
		X2[is.nan(X2)] <- 1
		X2[X2==Inf] <- 10^20

		Q1[i] <- as.numeric(quantile(X1,alpha))
		Q2[i] <- as.numeric(quantile(X2,alpha))
	
	}

	id <- 1:n
	
	R1 <- data.frame(id,VR,Q1)
	R1 <- R1[order(VR),]

	R2 <- data.frame(id,TR,Q2)
	R2 <- R2[order(TR),]

	return(list(VRATIO=R1,TAU2RATIO=R2))

}


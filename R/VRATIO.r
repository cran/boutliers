VRATIO <- function(y, v, method="REML", data=NULL, B=2000, alpha=0.05, seed=123456){

	set.seed(seed)

	if(is.null(data)==FALSE){

		data <- data.frame(data)

		y <- data[, deparse(substitute(y))]
		v <- data[, deparse(substitute(v))]

	}

	reml1 <- rma(yi=y,vi=v,method=method)
	
	mu0 <- as.numeric(reml1$beta)
	v0 <- v + reml1$tau2

	V0 <- reml1$tau2
	V1 <- reml1$se^2

	n <- length(y)
	VR <- TR <- numeric(n)

	for(i in 1:n){
	
		y_i <- y[setdiff(1:n,i)]
		v_i <- v[setdiff(1:n,i)]
	
		reml_i <- rma(yi=y_i,vi=v_i,method=method)
		
		VR[i] <- reml_i$se^2 / V1
		TR[i] <- reml_i$tau2 / V0		
		
	}
                              
	VR.b <- TR.b <- matrix(numeric(n*B),B)

	for(i in 1:n){
	
		y_i <- y[setdiff(1:n,i)]
		v_i <- v[setdiff(1:n,i)]
	
		reml_i <- rma(yi=y_i,vi=v_i,method=method)
	
		for(b in 1:B){
	
			y.b <- rnorm(n, mean=reml_i$beta, sd=sqrt(v+reml_i$tau2))
		
			y.b_i <- y.b[setdiff(1:n,i)]
	
			reml.b_1 <- rma(yi=y.b,vi=v,method=method)
			reml.b_0 <- rma(yi=y.b_i,vi=v_i,method=method)
		
			VR.b[b,i] <- reml.b_0$se^2 / (reml.b_1$se^2)
			TR.b[b,i] <- reml.b_0$tau2 / reml.b_1$tau2

		}
	
		print1 <- paste0("The bootstrap for ",i,"th study is completed.")
		if(b%%100==0) message(print1)

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


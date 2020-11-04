LRT_FE <- function(y, v, B=2000, alpha=0.05){

	ML0 <- ML_FE(y,v)

	mu0 <- ML0$mu
	v0 <- v + ML0$V0
	mlike0 <- ML0$Loglikelihood		# loglikelihood under H0

	n <- length(y)
	LR <- P <- numeric(n)

	for(i in 1:n){
	
		yi <- y[c(i,setdiff(1:n,i))]
		vi <- v[c(i,setdiff(1:n,i))]
	
		mlike1 <- SML_FE(yi,vi)$Loglikelihood
                              
		LR0 <- -2*(mlike0 - mlike1)		# LRT statistic

		LR[i] <- LR0
		P[i] <- 1 - pchisq(LR0,df=1)

	}
	
	LR.b <- matrix(numeric(n*B),B)

	for(b in 1:B){
	
		y.b <- rnorm(n, mean=mu0, sd=sqrt(v0))
		
		mlike0.b <- ML_FE(y.b,v)$Loglikelihood
		
		for(i in 1:n){
	
			yi <- y.b[c(i,setdiff(1:n,i))]
			vi <- v[c(i,setdiff(1:n,i))]
	
			mlike1.b <- SML_FE(yi,vi)$Loglikelihood
                              
			LR0.b <- -2*(mlike0.b - mlike1.b)		# LRT statistic
			LR.b[b,i] <- LR0.b

		}
	
		print1 <- paste0("The ",b,"th bootstrap is completed.")
		if(b%%100==0) print(print1)

	}
		
	P <- Q <- numeric(n)	
		
	for(i in 1:n){
	
		X.b <- LR.b[,i]
		P[i] <- QT(X.b, LR[i])
		Q[i] <- as.numeric(quantile(X.b,1-alpha))
	
	}

	id <- 1:n
	R <- data.frame(id,LR,Q,P)
	R <- R[order(P),]

	return(R)

}


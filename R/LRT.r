LRT <- function(y, v, model="RE", data=NULL, B=2000, alpha=0.05, seed=123456){

	set.seed(seed)

	if(is.null(data)==FALSE){

		data <- data.frame(data)

		y <- data[, deparse(substitute(y))]
		v <- data[, deparse(substitute(v))]

	}

	if(model=="RE"){

		ML0 <- ML(y,v)

		mu0 <- ML0$mu
		v0 <- v + ML0$V0
		mlike0 <- ML0$Loglikelihood		# loglikelihood under H0

		n <- length(y)
		LR <- P <- numeric(n)

		for(i in 1:n){
	
			yi <- y[c(i,setdiff(1:n,i))]
			vi <- v[c(i,setdiff(1:n,i))]
	
			mlike1 <- SML(yi,vi)$Loglikelihood
                              
			LR0 <- -2*(mlike0 - mlike1)		# LRT statistic

			LR[i] <- LR0
			P[i] <- 1 - pchisq(LR0,df=1)

		}
	
		LR.b <- matrix(numeric(n*B),B)

		for(i in 1:n){
	
			y_i <- y[setdiff(1:n,i)]
			v_i <- v[setdiff(1:n,i)]

			vi <- v[c(i,setdiff(1:n,i))]

			MLi0 <- ML(y_i,v_i)

			for(b in 1:B){
	
				y.b <- rnorm(n, mean=MLi0$mu, sd=sqrt(v+MLi0$V0))
				y.b <- y.b[c(i,setdiff(1:n,i))]
		
				mlike0.b <- ML(y.b,vi)$Loglikelihood
				mlike1.b <- SML(y.b,vi)$Loglikelihood
                              
				LR0.b <- -2*(mlike0.b - mlike1.b)		# LRT statistic
				LR.b[b,i] <- LR0.b

			}
	
			print1 <- paste0("The bootstrap for ",i,"th study is completed.")
			if(b%%100==0) message(print1)
		
		}
		
		P <- Q <- numeric(n)	
		
		for(i in 1:n){
	
			X.b <- LR.b[,i]
			P[i] <- QT(X.b, LR[i])
			Q[i] <- as.numeric(quantile(X.b,(1-alpha)))
	
		}

		id <- 1:n
		R <- data.frame(id,LR,Q,P)
		R <- R[order(P),]

		return(R)

	}
	
	if(model=="FE"){

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

		for(i in 1:n){
	
			y_i <- y[setdiff(1:n,i)]
			v_i <- v[setdiff(1:n,i)]

			vi <- v[c(i,setdiff(1:n,i))]

			MLi0 <- ML_FE(y_i,v_i)

			for(b in 1:B){
	
				y.b <- rnorm(n, mean=MLi0$mu, sd=sqrt(v+MLi0$V0))
				y.b <- y.b[c(i,setdiff(1:n,i))]
			
				mlike0.b <- ML_FE(y.b,vi)$Loglikelihood
				mlike1.b <- SML_FE(y.b,vi)$Loglikelihood
                              
				LR0.b <- -2*(mlike0.b - mlike1.b)		# LRT statistic
				LR.b[b,i] <- LR0.b

			}
	
			print1 <- paste0("The bootstrap for ",i,"th study is completed.")
			if(b%%100==0) message(print1)
		
		}
		
		P <- Q <- numeric(n)	
		
		for(i in 1:n){
	
			X.b <- LR.b[,i]
			P[i] <- QT(X.b, LR[i])
			Q[i] <- as.numeric(quantile(X.b,(1-alpha)))
	
		}

		id <- 1:n
		R <- data.frame(id,LR,Q,P)
		R <- R[order(P),]

		return(R)

	}

}


ML <- function(y,v,maxitr=200){

  N <- length(y)

  mu <- 0.1	# initial values
  V0 <- 0.1

  Qc0 <- c(mu,V0)

  LL1 <- function(V){

    ll1 <- 0

    for(i in 1:N){

      yi <- y[i]
      vi <- v[i]

      A1 <- 0.5*log(2*pi*(vi+V))
      A2 <- (yi-mu)^2/(2*(vi+V))

      ll1 <- ll1 + A1 + A2

    }

    return(ll1)

  }

  for(itr in 1:maxitr){

    wi <- (v + V0)^-1

    mu <- sum(wi * y)/sum(wi)
    V0 <- optimize(LL1, lower = 0, upper = 500)$minimum

    Qc <- c(mu,V0)

    rb <- abs(Qc - Qc0)/abs(Qc0); rb[is.nan(rb)] <- 0
    if(max(rb) < 10^-4) break

    Qc0 <- Qc

  }

  LL <- -LL1(V0)

  R1 <- list("mu"=mu,"V0"=V0,"Loglikelihood"=LL)

  return(R1)

}

ML_FE <- function(y,v,maxitr=200){

	N <- length(y)

	V0 <- 0

	LL1 <- function(V){

		ll1 <- 0
		
		for(i in 1:N){

			yi <- y[i]
			vi <- v[i]

			A1 <- 0.5*log(2*pi*(vi+V))
			A2 <- (yi-mu)^2/(2*(vi+V))
		
			ll1 <- ll1 + A1 + A2

		}

		return(ll1)
	
	}

	wi <- (v + V0)^-1
		
	mu <- sum(wi * y)/sum(wi)

	LL <- -LL1(V0)
	
	R1 <- list("mu"=mu,"V0"=V0,"Loglikelihood"=LL)
	
	return(R1)
	
}

QT <- function(x,x0){

 x1 <- sort(c(x,x0))
 w1 <- which(x1==as.numeric(x0))
 qt <- 1 - w1/(length(x)+1)
 return(qt)
 
}

REML <- function(y,v,maxitr=200){

	N <- length(y)

	mu <- 0.1	# initial values
	V0 <- 0.1
	
	Qc0 <- c(mu,V0)
	
	LL1 <- function(V){

		ll1 <- 0
		
		for(i in 1:N){

			yi <- y[i]
			vi <- v[i]

			A1 <- log(vi+V)
			A2 <- (yi-mu)^2/(vi+V)
		
			ll1 <- ll1 + A1 + A2

		}

		A3 <- log(sum((v+V)^-1))
		
		ll1 <- ll1 + A3

		return(ll1)
	
	}

	for(itr in 1:maxitr){
	
		wi <- (v + V0)^-1
		
		mu <- sum(wi * y)/sum(wi)
		V0 <- optimize(LL1, lower = 0, upper = 500)$minimum

		Qc <- c(mu,V0)

		rb <- abs(Qc - Qc0)/abs(Qc0); rb[is.nan(rb)] <- 0
		if(max(rb) < 10^-4) break
		
		Qc0 <- Qc
		
	}
	
	V1 <- sum((v + V0)^-1)^-1

	R1 <- list("mu"=mu,"V0"=V0,"V1"=V1)
	
	return(R1)
	
}

SML <- function(y,v,maxitr=200){

	N <- length(y)

	mu <- 0.1	# initial values
	beta <- 0.1
	V0 <- 0.1
	
	Qc0 <- c(mu,beta,V0)
	
	LL1 <- function(V){

		ll1 <- 0
		
		yi <- y[1]
		vi <- v[1]

		A1 <- 0.5*log(2*pi*(vi+V))
		A2 <- (yi-mu-beta)^2/(2*(vi+V))
		
		ll1 <- ll1 + A1 + A2

		for(i in 2:N){

			yi <- y[i]
			vi <- v[i]

			A1 <- 0.5*log(2*pi*(vi+V))
			A2 <- (yi-mu)^2/(2*(vi+V))
		
			ll1 <- ll1 + A1 + A2

		}

		return(ll1)
	
	}

	LL2 <- function(mu1){

		ll1 <- 0
		
		yi <- y[1]
		vi <- v[1]

		A1 <- 0.5*log(2*pi*(vi+V0))
		A2 <- (yi-mu1-beta)^2/(2*(vi+V0))
		
		ll1 <- ll1 + A1 + A2

		for(i in 2:N){

			yi <- y[i]
			vi <- v[i]

			A1 <- 0.5*log(2*pi*(vi+V0))
			A2 <- (yi-mu1)^2/(2*(vi+V0))
		
			ll1 <- ll1 + A1 + A2

		}

		return(ll1)
	
	}

	LL3 <- function(beta1){

		ll1 <- 0
		
		yi <- y[1]
		vi <- v[1]

		A1 <- 0.5*log(2*pi*(vi+V0))
		A2 <- (yi-mu-beta1)^2/(2*(vi+V0))
		
		ll1 <- ll1 + A1 + A2

		for(i in 2:N){

			yi <- y[i]
			vi <- v[i]

			A1 <- 0.5*log(2*pi*(vi+V0))
			A2 <- (yi-mu)^2/(2*(vi+V0))
		
			ll1 <- ll1 + A1 + A2

		}

		return(ll1)
	
	}

	for(itr in 1:maxitr){
	
		V0 <- optimize(LL1, lower = 0, upper = 500)$minimum
		mu <- optimize(LL2, lower = -500, upper = 500)$minimum
		beta <- optimize(LL3, lower = -500, upper = 500)$minimum
		
		Qc <- c(mu,beta,V0)

		rb <- abs(Qc - Qc0)/abs(Qc0); rb[is.nan(rb)] <- 0
		if(max(rb) < 10^-4) break
		
		Qc0 <- Qc
		
	}

	LL <- -LL1(V0)
	
	R1 <- list("mu"=mu,"beta"=beta,"V0"=V0,"Loglikelihood"=LL)
	
	return(R1)
		
}

SML_FE <- function(y,v,maxitr=200){

	N <- length(y)

	mu <- 0.1	# initial values
	beta <- 0.1
	V0 <- 0
	
	Qc0 <- c(mu,beta)
	
	LL1 <- function(V){

		ll1 <- 0
		
		yi <- y[1]
		vi <- v[1]

		A1 <- 0.5*log(2*pi*(vi+V))
		A2 <- (yi-mu-beta)^2/(2*(vi+V))
		
		ll1 <- ll1 + A1 + A2

		for(i in 2:N){

			yi <- y[i]
			vi <- v[i]

			A1 <- 0.5*log(2*pi*(vi+V))
			A2 <- (yi-mu)^2/(2*(vi+V))
		
			ll1 <- ll1 + A1 + A2

		}

		return(ll1)
	
	}

	LL2 <- function(mu1){

		ll1 <- 0
		
		yi <- y[1]
		vi <- v[1]

		A1 <- 0.5*log(2*pi*(vi+V0))
		A2 <- (yi-mu1-beta)^2/(2*(vi+V0))
		
		ll1 <- ll1 + A1 + A2

		for(i in 2:N){

			yi <- y[i]
			vi <- v[i]

			A1 <- 0.5*log(2*pi*(vi+V0))
			A2 <- (yi-mu1)^2/(2*(vi+V0))
		
			ll1 <- ll1 + A1 + A2

		}

		return(ll1)
	
	}

	LL3 <- function(beta1){

		ll1 <- 0
		
		yi <- y[1]
		vi <- v[1]

		A1 <- 0.5*log(2*pi*(vi+V0))
		A2 <- (yi-mu-beta1)^2/(2*(vi+V0))
		
		ll1 <- ll1 + A1 + A2

		for(i in 2:N){

			yi <- y[i]
			vi <- v[i]

			A1 <- 0.5*log(2*pi*(vi+V0))
			A2 <- (yi-mu)^2/(2*(vi+V0))
		
			ll1 <- ll1 + A1 + A2

		}

		return(ll1)
	
	}

	for(itr in 1:maxitr){
	
		mu <- optimize(LL2, lower = -500, upper = 500)$minimum
		beta <- optimize(LL3, lower = -500, upper = 500)$minimum
		
		Qc <- c(mu,beta)

		rb <- abs(Qc - Qc0)/abs(Qc0); rb[is.nan(rb)] <- 0
		if(max(rb) < 10^-4) break
		
		Qc0 <- Qc
		
	}

	LL <- -LL1(V0)
	
	R1 <- list("mu"=mu,"beta"=beta,"V0"=V0,"Loglikelihood"=LL)
	
	return(R1)
		
}




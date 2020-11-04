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



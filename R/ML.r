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

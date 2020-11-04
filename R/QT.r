QT <- function(x,x0){

 x1 <- sort(c(x,x0))
 w1 <- which(x1==as.numeric(x0))
 qt <- 1 - w1/(length(x)+1)
 return(qt)
 
}


convert_mean <- function(n1, m1, s1, n2, m2, s2, pooled = FALSE, type=c("MD","SMD")) {

  # initial check
  lstt <- c("MD", "SMD")
  type <- match.arg(type)
  
  # initial check
  util_check_nonneg(n1)
  #util_check_nonneg(m1)
  util_check_nonneg(s1)
  util_check_nonneg(n2)
  #util_check_nonneg(m2)
  util_check_nonneg(s2)
  util_check_num(pooled)

  if (length(n1) != length(m1) || length(m1) != length(s1) || length(s1) != length(n2) ||
      length(n2) != length(m2) || length(m2) != length(s2)) {
    stop("'n1', 'm1', 's1', 'n2', 'm2', and 's2' should have the same length.")
  }else if (!is.element(type, lstt)) {
    stop("Unknown 'type' specified.")
  }

  if(type=="MD"){
   res <- NULL
   res$y <- m1 - m2
   if (pooled) {
     res$v <- ((n1 - 1)*s2^2 + (n2 - 1)*s2^2)/(n1 + n2 - 2)*(1/n1 + 1/n2)
   } else {
     res$v <- s1^2/n1 + s2^2/n2
   }
  }
  
  if(type=="SMD"){
   res <- NULL
   s0 <- sqrt(((n1 - 1)*s1^2 + (n2 - 1)*s2^2)/(n1 + n2 - 2))
   res$y <- (m1 - m2)/s0
   res$v <- 1/n1 + 1/n2
  }

  res <- data.frame(res)

  return(res)

}

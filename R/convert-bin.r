convert_bin <- function(m1, n1, m2, n2, type = c("logOR", "logRR", "RD")) {

  # initial check
  lstt <- c("logOR", "logRR", "RD")
  type <- match.arg(type)

  util_check_nonneg(m1)
  util_check_nonneg(n1)
  util_check_nonneg(m2)
  util_check_nonneg(n2)

  if (length(m1) != length(n1) || length(n1) != length(m2) || length(m2) != length(n2)) {
    stop("'m1', 'n1', 'm2', and 'n2' should have the same length.")
  } else if (!is.element(type, lstt)) {
    stop("Unknown 'type' specified.")
  }

  res <- NULL
  if (type == "logOR") {
    res$y <- log((m1 + 0.5)*(n2 - m2 + 0.5)/(n1 - m1 + 0.5)/(m2 + 0.5))
    res$v <- 1.0/(m1 + 0.5) + 1.0/(n1 - m1 + 0.5) + 1.0/(m2 + 0.5) + 1.0/(n2 - m2 + 0.5)
  } else if (type == "logRR") {
    res$y <- log((m1 + 0.5)*(n2 + 0.5)/(n1 + 0.5)/(m2 + 0.5))
    res$v <- 1.0/(m1 + 0.5) - 1.0/(n1 + 0.5) + 1.0/(m2 + 0.5) - 1.0/(n2 + 0.5)
  } else if (type == "RD") {
    res$y <- m1/n1 - m2/n2
    res$v <- ((m1 + 0.0625)/(n1 + 0.125))*((n1 - m1 + 0.0625)/(n1 + 0.125))/n1 +
        ((m2 + 0.0625)/(n2 + 0.125))*((n2 - m2 + 0.0625)/(n2 + 0.125))/n2
  }

  res <- data.frame(res)

  return(res)

}

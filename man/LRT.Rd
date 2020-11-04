\name{LRT}
\alias{LRT}
\title{Likelihood ratio test using a mean-shifted model}
\usage{
LRT(y, v, B=2000, alpha=0.05)
}
\arguments{
\item{y}{A vector of the outcome measure estimates (e.g., MD, SMD, log OR, log RR, RD)}

\item{v}{A vector of the variance estimate of \code{y}}

\item{B}{The number of bootstrap resampling (defalt: 2000)}

\item{alpha}{The significance level (default: 0.05)}
}
\value{
Results of the likelihood ratio tests involving bootstrap p-values. The outputs are ordered by the p-values.
\itemize{
\item \code{id}: ID of the study.
\item \code{LR}: The likelihood ratio statistic for based on the mean-shifted model.
\item \code{Q}: \code{1-alpha}th percentile for the bootstrap distribution of the likelihood ratio statistic.
\item \code{P}: The bootstrap p-value for the likelihood ratio statistic.
}
}
\description{
Implementing the likelihood ratio tests using the mean-shifted model for the DerSimonian-Laird-type random-effects model. The bootstrap p-values are provided.
}
\examples{
data(SMT)
attach(SMT)

dat <- convert_mean(n1, m1, s1, n2, m2, s2, type="MD")

LRT(dat$y, dat$v, B=100)
# For quick demonstration, B is set to 100.
# Practically, B should be set to >= 1000 (default is 2000).

detach(SMT)
}

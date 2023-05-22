\name{LRT_FE}
\alias{LRT_FE}
\title{Likelihood ratio test using a mean-shifted model by the fixed-effect model}
\usage{
LRT_FE(y, v, B=2000, alpha=0.05)
}
\arguments{
\item{y}{A vector of the outcome measure estimates (e.g., MD, SMD, log OR, log RR, RD)}

\item{v}{A vector of the variance estimate of \code{y}}

\item{B}{The number of bootstrap resampling (default: 2000)}

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
Implementing the likelihood ratio tests using the mean-shifted model for the fixed-effect model. The bootstrap p-values are provided.
}
\examples{
require(metafor)
data(SMT)

edat2 <- escalc(m1i=m1,sd1i=s1,n1i=n1,m2i=m2,sd2i=s2,n2i=n2,measure="MD",data=SMT)

LRT_FE(edat2$yi, edat2$vi, B=10)
# This is an example command for illustration. B should be >= 1000.
}

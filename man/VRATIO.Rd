\name{VRATIO}
\alias{VRATIO}
\title{Variance ratio influential statistics}
\usage{
VRATIO(y, v, B=2000, alpha=0.05)
}
\arguments{
\item{y}{A vector of the outcome measure estimates (e.g., MD, SMD, log OR, log RR, RD)}

\item{v}{A vector of the variance estimate of \code{y}}

\item{B}{The number of bootstrap resampling (default: 2000)}

\item{alpha}{The bootstrap percentile to be outputted (default: 0.05)}
}
\value{
The variance ratio influential statistics by leave-one-out analysis and their bootstrap percentiles. The outputs are ordered by the sizes of the variance ratio statistics.
\itemize{
\item \code{id}: ID of the study.
\item \code{VR}: The VRATIO statistic (relative change of the variance of the overall estimator) by leave-one-out analysis. 
\item \code{Q1}: \code{alpha}th percentile for the bootstrap distribution of the VRATIO statistic.
\item \code{TR}: The TAU2RATIO statistic (relative change of the heterogeneity variance) by leave-one-out analysis. 
\item \code{Q2}: \code{alpha}th percentile for the bootstrap distribution of the TAU2RATIO statistic.
}
}
\description{
Calculating the variance ratio influential statistics by leave-one-out analysis and the percentiles of their bootstrap distributions.
}
\examples{
require(metafor)
data(finasteride)

edat3 <- escalc(m1i=m1,sd1i=s1,n1i=n1,m2i=m0,sd2i=s0,n2i=n0,measure="MD",data=finasteride)

VRATIO(edat3$yi, edat3$vi, B=10)
# This is an example command for illustration. B should be >= 1000.
}

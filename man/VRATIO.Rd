\name{VRATIO}
\alias{VRATIO}
\title{Variance ratio influential statistics}
\usage{
VRATIO(y, v, B=2000, alpha=0.05)
}
\arguments{
\item{y}{A vector of the outcome measure estimates (e.g., MD, SMD, log OR, log RR, RD)}

\item{v}{A vector of the variance estimate of \code{y}}

\item{B}{The number of bootstrap resampling (defalt: 2000)}

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
data(finasteride)
attach(finasteride)

dat <- convert_mean(n1, m1, s1, n0, m0, s0, type="MD")
print(dat)

VRATIO(dat$y, dat$v, B=100)
# For quick demonstration, B is set to 100.
# Practically, B should be set to >= 1000 (default is 2000).

detach(finasteride)
}

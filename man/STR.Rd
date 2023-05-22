\name{STR}
\alias{STR}
\title{Studentized residuals by leave-one-out analysis}
\usage{
STR(y, v, B=2000, alpha=0.95)
}
\arguments{
\item{y}{A vector of the outcome measure estimates (e.g., MD, SMD, log OR, log RR, RD)}

\item{v}{A vector of the variance estimate of \code{y}}

\item{B}{The number of bootstrap resampling (default: 2000)}

\item{alpha}{The bootstrap percentiles to be outputted; 0.5(1-alpha)th and (1-0.5(1-alpha))th pecentiles. Default is 0.95; 2.5th and 97.5th percentiles are calculated.}
}
\value{
The studentized residuals by leave-one-out analysis. The outputs are ordered by the sizes of the studentized residuals.
\itemize{
\item \code{id}: ID of the study.
\item \code{psi}: The studentized residuals by leave-one-out analysis (studentized deleted residuals).
\item \code{Q1}: 0.5(1-alpha)th percentile for the bootstrap distribution of the studentized residual (default: 2.5th percentile).
\item \code{Q2}: 1-0.5(1-alpha)th percentile for the bootstrap distribution of the studentized residual (default: 97.5th percentile).
}
}
\description{
Calculating the studentized residuals by leave-one-out analysis (studentized deleted residuals) and the percentiles of their bootstrap distributions.
}
\examples{
require(metafor)
data(PPI)

edat1 <- escalc(ai=d1,n1i=n1,ci=d2,n2i=n2,measure="OR",data=PPI)

STR(edat1$yi,edat1$vi,B=10)   # This is an example command for illustration. B should be >= 1000.
}

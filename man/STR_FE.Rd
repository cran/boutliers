\name{STR_FE}
\alias{STR_FE}
\title{Studentized residuals by leave-one-out analysis for the fixed-effect model}
\usage{
STR_FE(y, v, B=2000, alpha=0.95)
}
\arguments{
\item{y}{A vector of the outcome measure estimates (e.g., MD, SMD, log OR, log RR, RD)}

\item{v}{A vector of the variance estimate of \code{y}}

\item{B}{The number of bootstrap resampling (defalt: 2000)}

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
Calculating the studentized residuals by leave-one-out analysis (studentized deleted residuals) for the fixed-effect model and the percentiles of their bootstrap distributions.
}
\examples{
data(PPI)
attach(PPI)

dat <- convert_bin(d1, n1, d2, n2, type = "logOR")
print(dat)

STR_FE(dat$y, dat$v, B=100)
# For quick demonstration, B is set to 100.
# Practically, B should be set to >= 1000 (default is 2000).

detach(PPI)
}

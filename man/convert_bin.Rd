\name{convert_bin}
\alias{convert_bin}
\title{Converting binary data to summary statistics}
\usage{
convert_bin(m1, n1, m2, n2, type = c("logOR", "logRR", "RD"))
}
\arguments{
\item{m1}{A vector of the number of successes in experimental group}

\item{n1}{A vector of the number of patients in experimental group}

\item{m2}{A vector of the number of successes in contorol group}

\item{n2}{A vector of the number of patients in contorol group}

\item{type}{the outcome measure for binary outcome data (default = "logOR").
\itemize{
\item \code{logOR}: logarithmic odds ratio,
                    which is defined by \eqn{=\log \frac{(m1+0.5)(n2-m2+0.5)}{(n1-m1+0.5)(m2+0.5)}}.
\item \code{logRR}: logarithmic relative risk,
                    which is defined by \eqn{=\log \frac{(m1+0.5)(n2+0.5)}{(n1+0.5)(m2+0.5)}}.
\item \code{RD}:    risk difference,
                    which is defined by \eqn{=\frac{m1}{n1}-\frac{m2}{n2}}.
}}
}
\value{
A \code{data.frame} of study data.
\itemize{
\item \code{y}: A numeric vector of the effect size estimates.
\item \code{v}: A numeric vector of the within studies variances.
}
}
\description{
Converting binary outcome data to the effect size estimates
and the within studies variances vector
}
\details{
This function implements methods for logarithmic
odds ratio, logarithmic relative risk, and risk difference
described in Hartung & Knapp (2001).
}
\examples{
data(PPI)
attach(PPI)

dat1 <- convert_bin(d1, n1, d2, n2, type = "logOR")
print(dat1)

dat2 <- convert_bin(d1, n1, d2, n2, type = "logRR")
print(dat2)

dat3 <- convert_bin(d1, n1, d2, n2, type = "RD")
print(dat3)
}
\references{
Hartung, J., and Knapp, G. (2001).
A refined method for the meta-analysis of controlled clinical trials
with binary outcome.
\emph{Stat Med.}
\strong{20}(24): 3875-3889.
\url{https://doi.org/10.1002/sim.1009}
}

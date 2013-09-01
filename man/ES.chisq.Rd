\name{ES.chisq}
\alias{ES.chisq}
\title{Compute effect size of chi-squared test}
\usage{
  ES.chisq(ct = NULL, chisq = NULL, p = NULL, n = NULL,
    df = NULL, mindf = NULL)
}
\arguments{
  \item{ct}{a m x n Contingency Table (matrix with m rows
  and n columes)}

  \item{chisq}{the value the chi-squared test statistic}

  \item{p}{p value for the chi-squared test}

  \item{df}{degree of freedom (e.g., df=(m-1)*(n-1))}

  \item{n}{total number of observations (sample size)}

  \item{mindf}{the degrees of freedom for the variable with
  the smaller number of levels, if m > n, mindf=n-1,
  otherwise, mindf=m-1}
}
\value{
  es(effect size), chisq, p, df, n
}
\description{
  Compute effect size of chi-squared test
}
\examples{
counts <- matrix(c(225,125,85,95),nrow=2,byrow=TRUE);ES.chisq(ct=counts)
}
\seealso{
  \code{\link{power.chisq}}
}
\keyword{chi-squared}
\keyword{degree}
\keyword{effect}
\keyword{freedom}
\keyword{of}
\keyword{power,}
\keyword{size,}
\keyword{study}
\keyword{test,}


\name{power.chisq}
\alias{power.chisq}
\title{Compute study power of chi-squared test}
\usage{
  power.chisq(es = NULL, df = NULL, n = NULL, power = NULL,
    sig.level = NULL)
}
\arguments{
  \item{es}{effect size. A numeric value or output of
  ES.chisq}

  \item{df}{degree of freedom}

  \item{n}{total number of observations (sample size)}

  \item{power}{power of study}

  \item{sig.level}{significance level}
}
\value{
  power, es, df, n, sig.level
}
\description{
  Compute study power of chi-squared test
}
\examples{
counts <- matrix(c(225,125,85,95),nrow=2,byrow=TRUE);
power.chisq(es=ES.chisq(ct=counts),sig.level=0.05)
}
\seealso{
  \code{\link{ES.chisq}}
}
\keyword{chi-squared}
\keyword{effect}
\keyword{level}
\keyword{power,}
\keyword{sample}
\keyword{significant}
\keyword{size,}
\keyword{study}
\keyword{test,}


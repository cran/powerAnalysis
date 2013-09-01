\name{power.plot.chisq}
\alias{power.plot.chisq}
\title{Sample size estimation of chi-squared test}
\usage{
  power.plot.chisq(es = NULL, power = NULL, df = NULL,
    sig.level = NULL, allele = FALSE, xlab = "Effect Size",
    ylab = "Sample Size", main = "Sample Size Estimation",
    grid = FALSE, ...)
}
\arguments{
  \item{es}{effect size.}

  \item{power}{power of study}

  \item{df}{degree of freedom}

  \item{sig.level}{significance level}

  \item{allele}{in genetic association study, whether test
  allele or genotype}

  \item{xlab}{a title for the x axis}

  \item{ylab}{a title for the y axis}

  \item{main}{an overall title for the plot}

  \item{grid}{add grid lines or not}

  \item{...}{Arguments to be passed to methods}
}
\value{
  power, es, df, n, sig.level
}
\description{
  Sample size estimation of chi-squared test
}
\examples{
es=seq(from=0.05,to=0.5,by=0.05);
power=seq(from=0.7,to=0.9,by=0.1);
power.plot.chisq(es=es,power=power,df=1,sig.level=0.05,allele=TRUE)
}
\seealso{
  \code{\link{power.chisq}}
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


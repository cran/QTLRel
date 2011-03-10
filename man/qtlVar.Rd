\name{qtlVar}
\alias{qtlVar}
\title{
QTL Variance
}
\description{
Estimate variance in a quantitative trait induced by QTL.
}
\usage{
qtlVar(lrt,prdat,simulation=FALSE,nsim=25)
}
\arguments{
  \item{lrt}{
a data frame (a, d, g, e, ...). For instance, its components should be corresponding to those in an object of \code{\link{scanOne}} when polygenic variation is accounted for.
}
  \item{prdat}{
a 3-D array that provides probabilities of genotypes "AA", "AB" and "BB". If \code{prDat} is an object of \code{\link{genoProb}}, then \code{prdat} can be \code{prDat$pr}. 
}
  \item{simulation}{
whether to use simulations to estimate the variance explained by QTL.
}
  \item{nsim}{
number of simulations to perform if \code{simulation} is TRUE.
}
}
\value{
A vector displaying the estimated variance at each loci.
}
\note{
Correlations among observations are ignored, and this function should be used with caution.
}

\seealso{
\code{\link{scanOne}} and \code{\link{genoProb}}
}


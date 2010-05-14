\name{blup}
\Rdversion{0.1}
\alias{blup}
\title{
Best Linear Unbiased Prediction
}
\description{
Estimate the best linear unbiased prediction (BLUP) for various effects in the model.
}
\usage{
blup(object)
}
\arguments{
  \item{object}{
an object from \code{\link{estVC}} or \code{\link{aicVC}}.
}
}
\value{
\item{fixed}{BLUP for fixed effects.}
\item{AA,DD,\dots}{BLUP for random (genetic) variance components "AA", "DD", \dots}
\item{EE}{BLUP for residual effect.}
}

\seealso{
\code{\link{estVC}} and \code{\link{aicVC}}.
}
\examples{
y<- rnorm(100)
x<- matrix(1,nrow=100,ncol=1)
v<- cov(matrix(rnorm(10^6),ncol=100))

o<- estVC(y, x, v = list(AA=v,DD=NULL,HH=NULL,AD=NULL,
   MH=NULL,EE=diag(100)))
b<- blup(o)
}


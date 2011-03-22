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
a data frame (a, d, ...), where 'a' and 'd' are respectively additive and dominance effects.
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
\examples{
data(miscEx)

gdtmp<- (gdat=="AA") + (gdat=="AB")*2 + (gdat=="BB")*3
   gdtmp<- replace(gdtmp,is.na(gdtmp),0)
prDat<- genoProb(gdat=gdtmp, gmap=genMap, step=Inf,
   gr=2, method="Haldane", verbose=TRUE)
y<- rnorm(20)
x<- matrix(1,nrow=20,ncol=1)
v<- cov(matrix(rnorm(500*20),ncol=20))
o<- estVC(y, x, v = list(AA=v,DD=NULL,HH=NULL,AD=NULL,
   MH=NULL,EE=diag(20)))

gcv<- v*o$par["AA"]
llk.hk<- scanOne(y=y, x=x, vc=gcv, prdat=prDat)

qef<- NULL
for(n in 1:length(llk.hk$par))
   qef<- rbind(qef,llk.hk$par[[n]][c("a","d")])
   qef<- as.data.frame(qef)
qtlVar(qef,prDat$pr)[1:3]
}


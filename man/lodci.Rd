\name{lodci}
\alias{lodci}
\title{
Estimate LOD Support Intervals
}
\description{
Estimate LOD support intervals.
}
\usage{
lodci(llk,cv=0,lod=1.5,drop=3)
}
\arguments{
  \item{llk}{
a data frame with components (chr, dist, y, ...), where "chr" is the chromosome on which the scanning locus is located, "dist" is the genetic or physical position of the scanning locus, and "y" is the test statistic.
}
  \item{cv}{
threshold. Reported support intervals cover at least one scanning locus where \code{llk$y > cv}.
}
  \item{lod}{
\code{lod} (1.5 by default) LOD support intervals are reported when \code{llk$y} is converted to LOD score.
}
  \item{drop}{
3 by default. See "details".
}
}
\details{
In case of multiple peaks on a chromosome, a peak has to satisfy: a) above the threshold \code{cv}; b) drops, e.g., 3 LOD on both sides except chromosome ends. So if two peaks close to each other but LOD between them doesn't drop, e.g., 3 LOD, only one of them is considered.
}
\value{
A data frame (chr, lower, upper), where "chr" represents the chromosome, "lower" the lower bound and "upper" the upper bound.
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

tmp<- data.frame(y=llk.hk$p,chr=llk.hk$chr,dist=llk.hk$dist)
lodci(tmp, cv=5, lod=1.5, drop=3)
}


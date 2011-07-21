\name{plotit}
\alias{plot}
\alias{plotit}
\alias{plot.scanOne}
\alias{plot.scanTwo}
\title{
Plotting
}
\description{
Plot mapping results.
}
\usage{
\method{plot}{scanOne}(x,...)

plotit(lrt,cv,bychr=FALSE,chr.labels=TRUE,type="p",lty=NULL,
   col=NULL,pch=NULL,cex=NULL,...)
}
\arguments{
  \item{x}{
object from \code{\link{scanOne}} or \code{\link{scanTwo}}.
}
  \item{lrt}{
a data frame with (chr, dist, y,...) or (chr, dist, y, group,...), where "chr" represents chromosome, "dist" position on the chromosome, "y" the test statistic.
}
  \item{cv}{
threshold to be drawn on the plot.
}
  \item{cex}{
see \code{\link{par}}.
}
  \item{bychr}{
a logical variable. If true, the plot will be displayed per chromosomes.
}
  \item{chr.labels}{
a logical variable. If true, the chromosome names will be drawn.
}
  \item{type,lty,col,pch}{
See \code{\link{plot.default}}.
}
  \item{\dots}{
other options passed to R plot function. To call \code{\link{plot}} to plot results of \code{\link{scanOne}}, one may need to provide a genetic map \code{gmap} that should be data frame (snp, chr, dist, ...), where "snp" is the SNP (marker) name, "chr" is the chromosome where the "snp" is, and "dist" is the genetic distance in centi-Morgan (cM) from the leftmost SNP (marker) on the chromosome.

}
}
\note{
A genetic map 'gmap' may be needed to plot an object of \code{\link{scanOne}} or \code{\link{scanTwo}}. The color option may not give what is expected.
}
\examples{
data(miscEx)

\dontrun{
# impute missing genotypes
gdat.imp<- genoImpute(gdat, gmap=genMap, step=Inf,
   gr=2, na.str=NA)
y<- rnorm(20)
x<- matrix(1,nrow=20,ncol=1)
v<- cov(matrix(rnorm(500*20),ncol=20))
# estimate variance components
o<- estVC(y, x, v = list(AA=v,DD=NULL,HH=NULL,AD=NULL,
   MH=NULL,EE=diag(20)))

# genome scan
gcv<- v*o$par["AA"]
llk<- scanOne(y=y, x=x, vc=gcv, gdat=gdat.imp[,1:20])

# plotting
plot(llk,gmap=genMap) # gmap is needed

# plotting in another way
idx<- match(colnames(gdat.imp[,1:20]),genMap$snp)
tmp<- data.frame(chr=genMap$chr[idx],dist=genMap$dist[idx],y=llk$p)
plotit(tmp, main="Mapping Plot", xlab="Chromosome", ylab="LRT",
   col=as.integer(tmp$ch)\%\%2+2,type="b")
}
}


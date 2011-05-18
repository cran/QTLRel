\name{scanOne}
\alias{scanOne}
\title{
Genome Scan for QTL
}
\description{
Evaluate log-likelihood ratio test statistics or P-values at scanning loci along the genome.
}
\usage{
scanOne(y, x, gdat, prdat=NULL, vc=NULL, intcovar=NULL,
   numGeno=FALSE, test=c("None","F","Chisq"),
   minorGenoFreq=0, rmv=TRUE)
}
\arguments{
  \item{y}{
a numeric vector or a numeric matrix of one column (representing a phenotype for instance).
}
  \item{x}{
a data frame or matrix, representing covariates if not missing.
}
  \item{gdat}{
genotype data. Should be a matrix or a data frame, with each row representing an observation and each column a marker locus. The column names should be marker names. Optional if an object \code{prdat} from \code{\link{genoProb}} is used as an argument.
}
  \item{prdat}{
an object from \code{\link{genoProb}}.
}
  \item{vc}{
an object from \code{\link{estVC}} or \code{\link{aicVC}}, or an estimated variance-covariance matrix induced by relatedness and environment.
}
  \item{intcovar}{
Covariates that interact with QTL.
}
  \item{numGeno}{
whether to treat numeric coding of genotypes as numeric.
}
  \item{test}{
"None", "F" or "Chisq".
}
  \item{minorGenoFreq}{
specify the minimum tolerable minor genotype frequency at a scanning locus if \code{gdat} is used.
}
  \item{rmv}{
a logical variable. If true, then the scanning locus will be skipped if the minor genotype frequency at the locus is smaller than \code{minorGenoFreq}. Otherwise, the scanning process will stop and return with NULL.
}
}
\details{
The test at a scanning locus under the assumption of no QTL effect versus the assumption of QTL effect is performed by conditioning on the estimated polygenic genetic variance-covariance matrix. Normality is assumed for the random effects.
}
\value{
A list with the at least following components:
   \item{p}{P-value at the snp (marker) if \code{test} is "F" or "Chisq", or the log-likelihood ratio statistic at the SNP (marker) if "test" is "None"}
   \item{parameters}{estimated parameters at all scanning loci, including additive effect \code{a} and dominance effect \code{d} if \code{prdat} is not NULL}
}

\seealso{
\code{\link{genoImpute}} and  \code{\link{genoProb}}.
}
\examples{
data(miscEx)

gdat.imp<- genoImpute(gdat, gmap=genMap, step=Inf,
   gr=2, na.str=NA)
y<- rnorm(20)
x<- sample(c("Male","Female"),20,replace=TRUE)
v<- cov(matrix(rnorm(500*20),ncol=20))
o<- estVC(y, x, v = list(AA=v,DD=NULL,HH=NULL,AD=NULL,
   MH=NULL,EE=diag(20)))

pv<- scanOne(y=y, x=x, gdat=gdat.imp[,1:20], vc=o)
pv
plot(pv,gmap=genMap)

gdtmp<- (gdat=="AA") + (gdat=="AB")*2 + (gdat=="BB")*3
   gdtmp<- replace(gdtmp,is.na(gdtmp),0)
prDat<- genoProb(gdat=gdtmp, gmap=genMap, step=Inf,
   gr=2, method="Haldane", verbose=TRUE)

# Haley-Knott method
pv.hk<- scanOne(y=y, x=x, prdat=prDat, vc=o)
pv.hk
plot(pv.hk)
}


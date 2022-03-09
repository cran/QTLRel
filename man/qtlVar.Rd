\name{qtlVar}
\alias{qtlVar}
\title{
QTL Variance
}
\description{
Estimate variance in a quantitative trait induced by QTL.
}
\usage{
qtlVar(lrt, prdat, simulation = FALSE, nsim = 25)
}
\arguments{
  \item{lrt}{
A data frame (a, d, ...), where 'a' and 'd' are respectively additive and dominance effects.
}
  \item{prdat}{
A 3-D array that provides probabilities of genotypes "AA", "AB" and "BB". If \code{prDat} is an object of \code{\link{genoProb}}, then \code{prdat} can be \code{prDat$pr}. 
}
  \item{simulation}{
Whether to use simulations to estimate the variance explained by QTL.
}
  \item{nsim}{
Number of simulations to perform if \code{simulation} is TRUE.
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

\dontrun{
# impute missing genotypes
pheno<- pdatF8[!is.na(pdatF8$bwt) & !is.na(pdatF8$sex),]
ii<- match(rownames(pheno), rownames(gdatF8))
geno<- gdatF8[ii,]
ii<- match(rownames(pheno), rownames(gmF8$AA))
v<- list(A=gmF8$AA[ii,ii], D=gmF8$DD[ii,ii])

gdtmp<- geno
   gdtmp<- replace(gdtmp,is.na(gdtmp),0)
# rung 'genoProb'
prDat<- genoProb(gdat=gdtmp, gmap=gmapF8,
   gr=8, method="Haldane", msg=TRUE)
# estimate variance components
o<- estVC(y=pheno$bwt, x=pheno$sex, v=v)

# genome scan
pv.hk<- scanOne(y=pheno$bwt, x=pheno$sex, prdat=prDat, vc=o)

# run 'qtlVar'
qef<- pv.hk$par[,c("a","d")]
   qef<- as.data.frame(qef)
qv<- qtlVar(qef,prDat$pr)
}
}


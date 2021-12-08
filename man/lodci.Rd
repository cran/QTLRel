\name{lodci}
\alias{lodci}
\title{
Estimate LOD Support Intervals
}
\description{
Estimate LOD support intervals.
}
\usage{
lodci(llk, cv = 0, lod = 1.5, drop = 3)
}
\arguments{
  \item{llk}{
A data frame with components (chr, dist, y, ...), where "chr" is the chromosome on which the scanning locus is located, "dist" is the genetic or physical position of the scanning locus, and "y" is the test statistic.
}
  \item{cv}{
Threshold. Reported support intervals cover at least one scanning locus where \code{llk$y > cv}.
}
  \item{lod}{
LOD (specified by \code{lod}, which is 1.5 by default) support intervals to be reported when \code{llk$y} is converted to LOD score.
}
  \item{drop}{
3 by default. See "details".
}
}
\details{
In case of multiple peaks on a chromosome, a peak has to satisfy: a) above the threshold \code{cv}; b) drops, e.g., 3 LOD on both sides except chromosome ends. So if two peaks close to each other but LOD between them doesn't drop, e.g., 3 LOD, only one of them is considered.
}
\value{
A data frame with the following components:
\item{chr}{The chromosome}
\item{lower}{The lower bound}
\item{upper}{The upper bound}
\item{index}{Indicates which scanning loci}
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
# run 'genoProb'
prDat<- genoProb(gdat=gdtmp, gmap=gmapF8,
   gr=8, method="Haldane", msg=TRUE)
# estimate variance components
o<- estVC(y=pheno$bwt, x=pheno$sex, v=v)

# genome scan
llk.hk<- scanOne(y=pheno$bwt, x=pheno$sex, vc=o, prdat=prDat)

# extract LOD support intervals
tmp<- data.frame(y=llk.hk$LRT, chr=llk.hk$chr, dist=llk.hk$dist)
lodci(tmp, cv=10, lod=1.5, drop=3)
}
}


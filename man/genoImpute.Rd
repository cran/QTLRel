\name{genoImpute}
\Rdversion{0.1}
\alias{genoImpute}
\title{
Impute Genotypic Data
}
\description{
Impute missing genotypic data in advance intercross lines (AIL).
}
\usage{
genoImpute(gdat, gmap, prd=NULL, step=Inf, gr=2, pos=NULL,
   method=c("Haldane","Kosambi"), na.str="NA", verbose=FALSE)
}
\arguments{
  \item{gdat}{
genotype data. Should be a matrix or a data frame, with each row representing an observation and each column a marker locus. The column names should be marker names. Optional if an object \code{prd} from \code{\link{genoProb}} is used as an argument.
}
  \item{gmap}{
a genetic map. Should be data frame (snp, chr, dist,...), where "snp" is the SNP (marker) name, "chr" is the chromosome where the "snp" is, and "dist" is the genetic distance in centi-Morgan (cM) from the leftmost SNP (marker) on the chromosome.
}
  \item{prd}{
an object from \code{\link{genoProb}} if not \code{NULL}. See "details" for more information.
}
  \item{step}{
the maximum distance (in cM) between two adjacent loci for which the probabilities are calculated. The distance corresponds to the "cumulative" recombination rate at \code{gr}-th generation.
}
  \item{gr}{
the generation under consideration.
}
  \item{pos}{
data frame (chr, dist, snp, ...). If given, \code{step} will be ignored.
}
  \item{method}{
whether "Haldane" or "Kosambi" mapping function should be used.
}
  \item{na.str}{
string for missing values.
}
  \item{verbose}{
a logical variable. If TRUE, certain information will be printed out during calculation.
}
}
\details{
The missing genotypic value is randomly assigned with a probability conditional on the genotypes of the flanking SNPs (makers). 

An object, \code{prd}, from \code{\link{genoProb}} alone can be used for the purpose of imputation. Then, the output (especially the putative loci) will be determined by \code{prd}. Optionally, it can be used together with \code{gdat} so that missing values in \code{gdat} will be imputed if possible, depending on whether loci in the columns of \code{gdat} can be identified in the third dimension of \code{prd}; this won't change the original genotypic data. See examples.
}
\value{
A matrix with the number of rows being the same as \code{gdat} and with the number of columns depending on the SNP set in both \code{gdat} and \code{gmap} and the \code{step} length.
}
\note{
Currently only suitable for advanced intercross lines.
}

\seealso{
\code{\link{genoProb}}
}
\examples{
data(miscEx)

sum(is.na(gdat))
gdat[11:13,21:25]

gdtmp<- (gdat=="AA") + (gdat=="AB")*2 + (gdat=="BB")*3
   gdtmp<- replace(gdtmp,is.na(gdtmp),0)
prDat<- genoProb(gdat=gdtmp, gmap=genMap, step=Inf,
   gr=2, method="Haldane", verbose=TRUE)
prDat$pr[11:13,,21:25]

tmp<- genoImpute(prd=prDat)
dim(gdat)
sum(is.na(tmp))
tmp[11:13,21:25]

tmp<- genoImpute(gdat[1:15,],prd=prDat)
dim(gdat)
sum(is.na(tmp))
tmp[11:13,21:25]

tmp<- genoImpute(gdat, gmap=genMap, step=Inf,
   gr=2, na.str=NA)
dim(gdat)
sum(is.na(tmp))
tmp[11:13,21:25]

tmp<- genoImpute(gdat[1:15,], gmap=genMap, step=Inf,
   gr=2, na.str=NA, verbose=TRUE)
dim(gdat)
sum(is.na(tmp))
tmp[11:13,21:25]
}


\name{genoProb}
\alias{genoProb}
\title{
Probability of a Genotype.
}
\description{
Calculate the probability of a genotype at a locus conditional on the genotypes of its flanking markers in advance intercross lines (AIL).
}
\usage{
genoProb(gdat, gmap, step, gr = 2, pos = NULL, method=c("Haldane", "Kosambi"),
   msg = FALSE)
}
\arguments{
  \item{gdat}{
Genotype data. Should be a matrix or a data frame, with each row representing an observation and each column a marker locus. The column names should be marker names. Each entry should be 1, 2, 3 or 0, corresponding to "AA", "AB", "BB" or missing genotype.
}
  \item{gmap}{
A genetic map. Should be data frame (snp, chr, dist,...), where "snp" is the SNP (marker) name, "chr" is the chromosome where the "snp" is, and "dist" is the genetic distance in centi-Morgan (cM) from the left of the chromosome.
}
  \item{step}{
Optional. If specified, it is the maximum "cumulative" genetic distance (in cM) between two adjacent loci for which the probabilities are calculated. The genetic distance corresponds to the "cumulative" recombination rate at \code{gr}-th generation.
}
  \item{gr}{
The generation under consideration.
}
  \item{pos}{
Data frame (chr, dist, snp, ...). If given, \code{step} will be ignored.
}
  \item{method}{
Whether "Haldane" or "Kosambi" mapping function should be used.
}
  \item{msg}{
A logical variable. If TRUE, certain information will be printed out during calculation.
}
}
\details{
The "cumulative" genetic distance between any two adjacent loci for which probabilities are calculated is not larger than \code{step}. If \code{step} is missing or \code{step = Inf}, probabilities will only be calculated at loci in both the columns of \code{gdat} and the rows of \code{gmap}. If \code{step} is small, a large set of putative loci will be considered, including all loci defined by the columns of \code{gdat} and the rows of \code{gmap}. 
}
\value{
Probabilities for genotypes as well as genetic map information (snp,chr,dist)
\item{pr}{A 3-D array with the first dimension corresponding to that of \code{gdat}, the second to three genotype and the third to the putative loci. The probabilities will be -1 if not imputable, which happens when the genotype data is missing at all loci on the chromosome.}
}
\note{
Currently only suitable for advanced intercross lines.
}
\examples{
data(miscEx)

\dontrun{
# briefly look at genotype data
sum(is.na(gdatF8))
gdatF8[1:5,1:5]

gdtmp<- gdatF8
   gdtmp<- replace(gdtmp,is.na(gdtmp),0)
# In case an individual is not imputable, then
# one needs to assign genotypes manually
prDat<- genoProb(gdat=gdtmp, gmap=gmapF8, gr=8, method="Haldane", msg=TRUE)
prDat$pr[1:5,,1:5]
}
}


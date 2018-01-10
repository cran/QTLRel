\name{hapSim}
\alias{hapSim}
\title{
Generate Genotypic Data
}
\description{
Simulate gametic data from a pedigree.
}
\usage{
hapSim(ped, gmap, ids, hap, method = c("Haldane", "Kosambi"))
}
\arguments{
  \item{ped}{
A pedigree, which is a data frame (id, sex, sire, dam, ...). In "sex", male should be "M", "Male" or 1, and female should be "F", "Female" or 2 (other than 0 and 1). If given, "generation" can be numeric 0, 1, 2, ... or non-numeric "F0", "F1", "F2", ..., which should be in an increasing order. Note that 0 is reserved for missing values. If a sire/dam is an inbred founder, its ID should be tagged by character 'i' (e.g. 1i, 2i, etc.). See \code{\link{pedRecode}}.
}
  \item{gmap}{
A genetic map. Should be data frame (snp, chr, dist, ...), where "snp" is the SNP (marker) name, "chr" is the chromosome where the "snp" is, and "dist" is the genetic distance in centi-Morgan (cM) from the left of the chromosome. If \code{gmap} is missing but \code{hap} not, all but the first two columns of \code{hap} are ignored.
}
  \item{ids}{
Genotypic data are extracted only for individuals with IDs specified by \code{ids}. If missing, genotypic data are extracted for all individuals in the pedigree. If \code{ped} is an object of \code{\link{pedRecode}}, \code{ids} should be referred to "old" IDs.
}
  \item{hap}{
Founders' haplotype data if not missing. Rows correspond to founders as specified by row names, and columns correspond to loci in the genetic map \code{gmap} in the exact order. For an individual, the haplotype should be (f1 m1 f2 m2 ...) where fi is the allele from father at the i-th locus and mi is the allele from mother at the i-th locus. Elements should be non-negative integers that are not larger than 16384. If \code{hap} is not supplied, founders are assumed to be inbred.
}
  \item{method}{
Whether "Haldane" or "Kosambi" mapping function should be used. This will be ignored if the recombination rate \code{recRate} is a component of \code{gmap}.
}
}
\details{
The pedigree should be in the same format as an output of \code{\link{pedRecode}}. Founders mean those whose parents
have 0 or negative IDs after the pedigree is recoded by \code{\link{pedRecode}}.
}
\value{
A matrix giving haplotypes.
}
\seealso{
\code{\link{pedRecode}} for more information.
}
\examples{
data(miscEx)

\dontrun{
# prepare pedigree in desired format
pedR<- pedRecode(pedF8)
pedR[1:5,] # check to find out three founders
# fake founder haplotypes
hapDat<- rbind(rep(1:2,nrow(gmapF8)),rep(3:4,nrow(gmapF8)),rep(5:6,nrow(gmapF8)))
rownames(hapDat)<- c("32089","1","2")
# simulate hyplotypes for F8 individuals
hd<- hapSim(pedF8, gmapF8, ids=pedF8$id[pedF8$gen=="F8"], hap=hapDat)
dim(hd)
hd[1:5,1:10]
}
}


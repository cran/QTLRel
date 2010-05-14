\name{nullSim}
\Rdversion{0.1}
\alias{nullSim}
\title{
Simulate null distribution
}
\description{
Simulate distribution under the hypothesis of no QTL by permutation (of genotypic data) or gene dropping.
}
\usage{
nullSim(y,x,gdat,ped,gmap,ids,method=c("permutation","gene dropping"),
   vc=NULL,intcovar=NULL,minorGenoFreq=0.05,rmv=TRUE,nit=25,ntimes=10)
}
\arguments{
  \item{y}{
a numeric vector or a numeric matrix of one column (representing a phenotype for instance).
}
  \item{x}{
a data frame or matrix, representing covariates if not missing.
}
  \item{gdat}{
genotype data without missing values. Should be a matrix or a data frame, with each row representing an observation and each column a marker locus. Ignored in the case of gene dropping.
}
  \item{ped}{
a pedigree, which is a data frame (id, sex, sire, dam, ...). In "sex", male should be "M", "Male" or 1, and female should be "F", "Female" or 2 (other than 0 and 1). If given, "generation" can be numeric 0, 1, 2, ... or non-numeric "F0", "F1", "F2", ... Note that 0 is reserved for missing values. Ignored in the case of permutation.
}
  \item{gmap}{
a genetic map. Should be data frame (snp, chr, dist, ...), where "snp" is the SNP (marker) name, "chr" is the chromosome where the "snp" is, and "dist" is the genetic distance in centi-Morgan (cM) from the leftmost SNP (marker) on the chromosome. Ignored in the case of permutation.
}
  \item{ids}{
IDs of the individuals in the data, whose values are given by \code{y} in the order. Ignored in the case of permutation.
}
  \item{method}{
permutation or gene dropping.
}
  \item{vc}{
an object from \code{\link{estVC}} or \code{\link{aicVC}}, or an estimated variance-covariance matrix induced by relatedness. The scan will assume no polygenic variation if \code{vc} is NULL.
}
  \item{intcovar}{
Covariates that interact with QTL.
}
  \item{minorGenoFreq}{
specify the minimum tolerable minor genotype frequency at a scanning locus if \code{gdat} is used.
}
  \item{rmv}{
a logical variable. If true, then the scanning locus will be skipped if the minor genotype frequency at the locus is smaller than \code{minorGenoFreq}. Otherwise, the scanning process will stop and return with NULL.
}
  \item{nit}{
number of iterations to call \code{\link{optim}} for optimization. Ignored if \code{vc} is NULL.
}
  \item{ntimes}{
number of simulations.
}
}
\value{
A vector of numbers of length \code{ntimes}.
}
\note{
Not intended for interval mapping.
}


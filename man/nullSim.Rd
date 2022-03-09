\name{nullSim}
\alias{nullSim}
\title{
Simulate null distribution
}
\description{
Simulate the distribution of the test statistic by permutation (of genotypic data) or gene dropping.
}
\usage{
nullSim(y, x, gdat, prdat, ped, gmap, hap,
	method = c("permutation","gene dropping"), vc = NULL, intc = NULL,
	numGeno = FALSE, test = c("None","F","LRT"), minorGenoFreq = 0.05,
	rmv = TRUE, gr = 2, ntimes = 1000)
}
\arguments{
  \item{y}{
A numeric vector or a numeric matrix of one column (representing a phenotype for instance).
}
  \item{x}{
A data frame or matrix, representing covariates if not missing.
}
  \item{gdat}{
Genotype data without missing values. Should be a matrix or a data frame, with each row representing a sample and each column a marker locus. Ignored in the case of gene dropping.
}
  \item{prdat}{
An object from \code{\link{genoProb}}, or in the same form.
}
  \item{ped}{
A pedigree, which is a data frame (id, sex, father/sire, mother/dam, ...). In "sex", male should be "M", "Male" or 1, and female should be "F", "Female" or 2 (other than 0 and 1). If given, "generation" can be numeric 0, 1, 2, ... or non-numeric "F0", "F1", "F2", ... Note that 0 is reserved for missing values. Ignored in the case of permutation.
}
  \item{gmap}{
A genetic map. Should be data frame (snp, chr, dist, ...), where "snp" is the SNP (marker) name, "chr" is the chromosome where the "snp" is, and "dist" is the genetic distance in centi-Morgan (cM) from the left of the chromosome. Ignored in the case of permutation.
}
  \item{hap}{
Founders' haplotype data if not missing. Rows correspond to all founders, which should be in the first places in the pedigree \code{ped} and in the exact order, and columns correspond to loci in the genetic map \code{gmap} in the exact order. For a sample, the haplotype should be (f1 m1 f2 m2 ...) where fi is the allele from father at the i-th locus and mi is the allele from mother at the i-th locus. Elements should be non-negative integers that are not larger than 16384. If missing, two founders with alleles 1 and 2 are assumed.
}
  \item{method}{
Permutation or gene dropping.
}
  \item{vc}{
An object from \code{\link{estVC}} or \code{\link{aicVC}}, or an estimated variance-covariance matrix induced by relatedness. The scan will assume no polygenic variation if \code{vc} is NULL if any locci have a genotype frequency smaller than \code{minorGenoFreq}.
}
  \item{intc}{
Covariates that interact with QTL.
}
  \item{numGeno}{
Whether to treat numeric coding of genotypes as numeric. If true, \code{minorGenoFreq} will be ignored.
}
  \item{test}{
"None", "F" or "LRT". Note: the result will be on the scale of -log10(p-value) if the test is "F" or "LRT"; otherwise, the result will be the log-likelihood test statistic. Moerover, the result from each simulation is the maximum over all the SNPs/variants. Therefore, the user should make sure what is pertinent.
}
  \item{minorGenoFreq}{
Specify the minimum tolerable minor genotype frequency at a scanning locus if \code{gdat} is used.
}
  \item{rmv}{
A logical variable. If true, then the scanning locus will be skipped if the minor genotype frequency at the locus is smaller than \code{minorGenoFreq}. Otherwise, the scanning process will stop and return with NULL.
}
  \item{gr}{
The generation under consideration.
}
  \item{ntimes}{
Number of simulations.
}
}
\details{
Two methods considered here are permutation test and gene dropping test as described as follows.

Permutation test: Depending on the genome-scan, one can provide either \code{gdat} or \code{prdat} respectively corresponding to single-marker analysis or interval mapping. Then only arguments in \code{\link{scanOne}} are needed in addition to \code{method} and \code{ntimes}.

Gene dropping test: If \code{prdat} is provided, then \code{gdat} will be ignored. The procedure will first call \code{\link{genoSim}} to generate new genotype data and then call \code{\link{genoProb}} to generate data for Haley-Knott interval mapping. If \code{prdat} is not provided, then \code{gdat} should be provided. The procedure will generate new genotype data and scan the genome using these generated genotype data. Haldane mapping function is used to generate data.
}
\value{
A vector of length \code{ntimes}, the n-th element of which is maximum of the test statistics (LRT or -log10(p-value)) over the n-th genome scan.
}
\seealso{
\code{\link{genoSim}},  \code{\link{genoProb}} and \code{\link{scanOne}}.
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

gdatTmp<- genoImpute(geno, gmap=gmapF8,
   gr=8, na.str=NA)
# estimate variance components
o<- estVC(y=pheno$bwt, x=pheno$sex, v=v)

# scan marker loci & permutation
ex1<- nullSim(y=pheno$bwt, x=pheno$sex, gdat=gdatTmp,
	method="permutation", vc=o, ntimes=10)

# Haley-Knott method & permutation
gdtmp<- geno
   gdtmp<- replace(gdtmp,is.na(gdtmp),0)
prDat<- genoProb(gdat=gdtmp, gmap=gmapF8,
   gr=8, method="Haldane", msg=TRUE)
ex2<- nullSim(y=pheno$bwt, x=pheno$sex, prdat=prDat,
	method="permutation", vc=o, ntimes=10)

# remove samples whose father is troublesome "32089" 
#    before running gene dropping
# otherwise, "hap" data needs to be supplied

# scan marker loci & gene dropping
idx<- is.element(rownames(pdatF8), pedF8$id[pedF8$sire=="32089"])
pheno<- pdatF8[!is.na(pdatF8$bwt) & !is.na(pdatF8$sex) & !idx,]
ii<- match(rownames(pheno), rownames(gdatF8))
geno<- gdatF8[ii,]
ii<- match(rownames(pheno), rownames(gmF8$AA))
v<- list(A=gmF8$AA[ii,ii], D=gmF8$DD[ii,ii])

gdatTmp<- genoImpute(geno, gmap=gmapF8,
   gr=8, na.str=NA)
o<- estVC(y=pheno$bwt, x=pheno$sex, v=v)

ex3<- nullSim(y=pheno$bwt, x=pheno$sex, gdat=gdatTmp, ped=pedF8,
	gmap=gmapF8, method="gene", vc=o, ntimes=10)

# Haley-Knott method & gene dropping
gdtmp<- geno
   gdtmp<- replace(gdtmp,is.na(gdtmp),0)
prDat<- genoProb(gdat=gdtmp, gmap=gmapF8,
   gr=8, method="Haldane", msg=TRUE)
ex4<- nullSim(y=pheno$bwt, x=pheno$sex, prdat=prDat, ped=pedF8,
	gmap=gmapF8, method="gene", vc=o, gr=8, ntimes=10)
}
}


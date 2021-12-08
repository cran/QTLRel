\name{scanOne}
\alias{scanOne}
\title{
Genome Scan for QTL
}
\description{
Likelihood ratio tests or F tests at scanning loci over the genome.
}
\usage{
scanOne(y, x, gdat, prdat = NULL, vc = NULL, intc = NULL,
   numGeno = FALSE, test = c("None","F","LRT"),
   minorGenoFreq = 0, rmv = TRUE)
}
\arguments{
  \item{y}{
A numeric vector or a numeric matrix of one column, representing a phenotype.
}
  \item{x}{
A data frame or matrix, representing covariates if not missing.
}
  \item{gdat}{
Genotype data. It should be a matrix or a data frame, with each row being a sample and each column a locus. The column names should be marker names. It is ignored if an object \code{prdat} from \code{\link{genoProb}} is used as an argument.

If \code{gdat} is not numeric, there can be more than three genotypes and all scanning loci should have the same number of genotypes.
}
  \item{prdat}{
An object from \code{\link{genoProb}}, or in the same form. It should have a class "addEff" if allelic effects are assumed to be additive (see example below).
}
  \item{vc}{
An object from \code{\link{estVC}} or \code{\link{aicVC}}, or an estimated variance-covariance matrix induced by relatedness and environment.
}
  \item{intc}{
Covariates that interact with QTL.
}
  \item{numGeno}{
Whether to treat numeric coding of genotypes as numeric. If true, \code{minorGenoFreq} will be ignored.
}
  \item{test}{
"None", "F" or "LRT".
}
  \item{minorGenoFreq}{
Specify the minimum tolerable minor genotype frequency at a scanning locus if \code{gdat} is used.
}
  \item{rmv}{
A logical variable. If true, then the scanning locus will be skipped if the minor genotype frequency at the locus is smaller than \code{minorGenoFreq}. Otherwise, the scanning process will stop and return with NULL if any locci have a genotype frequency smaller than \code{minorGenoFreq}.
}
}
\details{
The test at a scanning locus under the null hypothesis of no QTL effect is performed by conditioning on the polygenic genetic variance-covariance. Normality is assumed for the random effects.

It is possible to extend the Haley-Knott approach to multiple-allelic cases under the assumption that allelic effects are all additive. Then, \code{prdat} should be provided and be of class "addEff".
}
\value{
A list with at least the following components:
   \item{F or LRT}{the F-test or likelihood ratio test (LRT) statistic at the SNP (marker) if \code{test} is "F" or otherwise}
   \item{pval}{P-value at the snp (marker) if \code{test} is "F" or "LRT"}
   \item{v}{Variation explained by the SNP (marker)}
   \item{parameters}{Estimated parameters at all scanning loci, including additive effect \code{a} and dominance effect \code{d} if \code{prdat} is not NULL}
}
\references{
Haley, C. S., and S. A. Knott (1992). A simple regression method for mapping quantitative trait loci in line crosses using flanking markers. Heredity 69: 315-324.
}

\seealso{
\code{\link{genoImpute}} and  \code{\link{genoProb}}.
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

# estimate variance components
o<- estVC(y=pheno$bwt, x=pheno$sex, v=v)

# impute missing genotypes
gdtmp<- genoImpute(geno, gmap=gmapF8, gr=8, na.str=NA, msg=FALSE)
# genome scan and plotting
lrt<- scanOne(y=pheno$bwt, x=pheno$sex, gdat=gdtmp, vc=o)
lrt
plot(lrt,gmap=gmapF8)

# Haley-Knott method
gdtmp<- geno; unique(unlist(gdtmp))
   gdtmp<- replace(gdtmp,is.na(gdtmp),0)
prDat<- genoProb(gdat=gdtmp, gmap=gmapF8, gr=8, method="Haldane", msg=TRUE)
pv.hk<- scanOne(y=pheno$bwt, intc=pheno$sex, prdat=prDat, vc=o, test="F")
pv.hk
plot(pv.hk, gmap=gmapF8)

# assume additive allelic effects
class(prDat)<- c(class(prDat), "addEff")
lrt.hk<- scanOne(y=pheno$bwt, intc=pheno$sex, prdat=prDat, vc=o)
lrt.hk
}
}


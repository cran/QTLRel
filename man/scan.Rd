\name{scanOne}
\Rdversion{0.1}
\alias{scanOne}
\title{
Genome Scan
}
\description{
Evaluate maximum likelihood ratio test statistics at scanning loci along the genome.
}
\usage{
scanOne(y, x, gdat, prdat=NULL, vc=NULL, intcovar=NULL,
   minorGenoFreq=0, rmv=TRUE, nit=25)
}
\arguments{
  \item{y}{
a numeric vector or a numeric matrix of one column (representing a phenotype for instance).
}
  \item{x}{
a data frame or matrix, representing covariates if not missing.
}
  \item{gdat}{
genotype data without missing values. Should be a matrix or a data frame, with each row representing an observation and each column a marker locus. The column names should be marker names. Ignored if an object \code{prdat} from \code{\link{genoProb}} is used as an argument and \code{vc} is not NULL.
}
  \item{prdat}{
an object from \code{\link{genoProb}}. Ignored if \code{vc} is NULL.
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
}
\details{
Maximum likelihood is evaluated at a scanning locus. If polygenic variation is assumed, the evaluation is approximated by conditioning on estimated polygenic genetic variance-covariance matrix.

Normality is assumed for the random effects. The optimization function \code{\link{optim}} is adopted in the above function to estimate the parameters and maximum likelihood. The maximum likelihood ratio test statistic is defined as -2\{log(H0)-log(Ha)\}, where log(H0) is the log-likelihood with the assumption of no QTL effect and log(Ha) is the log-likelihood with the assumption of QTL effect.

If \code{prdat} is specified, then "Haley-Knott" regression method will be implemented.
}
\value{
A list with the following components:
   \item{lr}{maximum likelihood ratio test statistic at the snp (marker)}
   \item{snp}{snp (marker) name}
   \item{parameters}{estimated parameters at all scanning loci}
}
\references{
Haley, C. S., and S. A. Knott (1992). A simple regression method for mapping quantitative trait loci in line crosses using flanking markers. Heredity 69: 315-324.
}

\seealso{
\code{\link{estVC}}, \code{\link{scanOne4p}}.
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

llk<- scanOne(y=y, x=x, gdat=gdat.imp[,1:20], vc=o)
llk
plot(llk,gmap=genMap)

gdtmp<- (gdat=="AA") + (gdat=="AB")*2 + (gdat=="BB")*3
   gdtmp<- replace(gdtmp,is.na(gdtmp),0)
prDat<- genoProb(gdat=gdtmp, gmap=genMap, step=Inf,
   gr=2, method="Haldane",verbose=TRUE)

# Haley-Knott method
llk.hk<- scanOne(y=y, prdat=prDat, vc=o, intcovar=x)
llk.hk
plot(llk.hk)

# ignore polygenic variation
llk.0<- scanOne(y=y, gdat=gdat.imp[,1:20],intcovar=x)
llk.0
}


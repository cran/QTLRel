\name{cic}
\alias{cic}
\title{
Calculate Jacquard condensed identity coefficients
}
\description{
Calculate Jacquard condensed identity coefficients from a pedigree.
}
\usage{
cic(ped,ids,inter,df=2,ask=FALSE,verbose=FALSE)
}
\arguments{
  \item{ped}{
a pedigree, which is a data frame (id, sire, dam, ...). If given, "generation" can be numeric 0, 1, 2, ... or non-numeric "F0", "F1", "F2", ... If "sex" is included, male should be "M", "Male" or 1, and female should be "F", "Female" or 2 (other than 0 and 1). Note that 0 is reserved for missing values.
}
  \item{ids}{
IDs of the individuals for which to calculate the Jacquard condensed identity coefficients. If missing, all individuals in the pedigree \code{ped} will be considered.
}
  \item{inter}{
if given, it should indicate intermediate generations.
}
  \item{df}{
if \code{inter} is missing, \code{df} is used to derive (optimal) \code{inter}. If \code{df = 0}, then there will no intermediate generations. If \code{df} is large, then all generations will be used as intermediate generations.
}
  \item{ask}{
if true, users will be asked whether to proceed.
}
  \item{verbose}{
if true, will print out some messages.
}
}
\details{
The coefficients will be calculated for individuals with IDs specified by \code{ids}. All individuals will be considered if \code{ids} is missing. This is not recommended if the total number of individuals in the pedigree is large. Instead, it is recommended that \code{ids} is specified for interested individuals only

\code{df} is a tuning parameter. It should not be 0 (or smaller than 1) if the pedigree is large in depth (many generations) but the number of individuals is not small; otherwise, it can take forever to finish. It should not be \code{Inf} (or a large number) if the number of individuals in certain intermediate generation is very large.
}
\value{
A matrix G with G[,j] being the j-th Jacquard identity coefficients.
}
\references{
Abney, M., M. S. McPeek, and C. Ober (2000). Estimation of variance components of quantitative traits in inbred populations. Am. J. Hum. Genet. 141, 629-650.
}
\note{
It is better to remove the working directory if the program is interrupted by external forces (e.g. killed by users).

It may require your operating system support "long long" integer type in C++. If you run this function in a windows system, make sure the working directory is under system volume C and you have the write privilege.

Warning: you may need to run this program on a 64-bit machine in case of seeing such a message!
}
\examples{
data(miscEx)

id<- sample(ped$id[300:500],20)

\dontrun{
oo<- cic(ped,id,df=Inf,verbose=TRUE)
o0<- cic(ped,id,df=0,verbose=TRUE)
o2<- cic(ped,id,ask=TRUE,verbose=TRUE)
c(sum(abs(oo-o0) >1e-7),sum(abs(o2-o0) >1e-7)) # same results
}
}

\name{genMatrix}
\alias{genMatrix}
\title{
Derive genetic matrices
}
\description{
Derive genetic matrices from Jacquard condensed identity coefficients or genotypic data.
}
\usage{
genMatrix(x)
}
\arguments{
  \item{x}{
An object of \code{\link{cic}} or \code{\link{ibs}}, or genotypic data in a matrix or a data frame with each row representing an individual and each column a marker locus and entry being "AA", "AB", "BB" (or 1, 2, 3) without missing genotypes.
}
}
\value{
\item{AA}{Additive genetic matrix.}
\item{DD}{Dominance genetic matrix.}
\item{AD,HH,MH}{Other three genetic matrices (see Abney et. al. 2000).}
\item{ib}{Inbreeding coefficients.}
}
\references{
Abney, M., M. S. McPeek, and C. Ober (2000). Estimation of variance components of quantitative traits in inbred populations. Am. J. Hum. Genet. 141, 629-650.
}
\seealso{
\code{\link{cic}}
}
\examples{
data(miscEx)

ids<- sample(pedF8$id[300:500],20)

\dontrun{
# get condensed identity coefficients
oo<- cic(pedF8, ids=ids, df=0)
ksp<- kinship(pedF8, ids=ids) # kinship coefficients only
# extract genetic matrices
gm<- genMatrix(oo)
sum((gm$AA-2*ksp)>1e-7) # same results
}
}

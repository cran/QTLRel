\name{genMatrix}
\alias{genMatrix}
\title{
Derive genetic matrices
}
\description{
Derive genetic matrices from Jacquard condensed identity coefficients.
}
\usage{
genMatrix(idcf)
}
\arguments{
  \item{idcf}{
\code{idcf} an object of \code{\link{cic}}.
}
}
\value{
\item{AA}{additive genetic matrix.}
\item{DD}{dominance genetic matrix.}
\item{AD,HH,MH}{other three genetic matrices (see Abney et. al. 2000).}
\item{ib}{inbreeding coefficients.}
}
\references{
Abney, M., M. S. McPeek, and C. Ober (2000). Estimation of variance components of quantitative traits in inbred populations. Am. J. Hum. Genet. 141, 629-650.
}
\seealso{
\code{\link{cic}}
}
\examples{
data(miscEx)

id<- sample(ped$id[300:500],20)
oo<- cic(ped,id,df=0)

ksp<- kinship(ped,ids=id)
y<- genMatrix(oo)
sum((y$AA-2*ksp)>1e-7) # same results
}

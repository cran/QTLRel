\name{blup}
\alias{blup}
\title{
Best Linear Unbiased Prediction
}
\description{
Estimate the best linear unbiased prediction (BLUP) for various effects in the model.
}
\usage{
blup(object)
}
\arguments{
  \item{object}{
An object from \code{\link{estVC}} or \code{\link{aicVC}}.
}
}
\value{
\item{fixed}{BLUP for fixed effects.}
\item{R, etc.}{BLUP for random effects.}
}

\seealso{
\code{\link{estVC}} and \code{\link{aicVC}}.
}
\examples{
data(miscEx)

\dontrun{
# only consider additive genetic variance component
pheno<- pdatF8[!is.na(pdatF8$bwt) & !is.na(pdatF8$sex),]
ii<- match(rownames(pheno), rownames(gmF8$AA))
v<- list(A=gmF8$AA[ii,ii],D=gmF8$DD[ii,ii])
vc<- estVC(y=pheno$bwt, x=pheno$sex, v=v)
b<- blup(vc)
}
}


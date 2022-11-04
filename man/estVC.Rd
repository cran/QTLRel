\name{estVC}
\alias{estVC}
\title{
Estimate Variance Component Parameters
}
\description{
Estimate model parameters for covariates, genetic variance components and residual effect.
}
\usage{
estVC(y, x, v = list(E=diag(length(y))), initpar, nit = 25,
   method = c("ML", "REML"), control = list(), hessian = FALSE)
}
\arguments{
  \item{y}{
A numeric vector or a numeric matrix of one column (representing a phenotype for instance).
}
  \item{x}{
A data frame or matrix, representing covariates if not missing.
}
  \item{v}{
A list of matrices representing variance components of interest. Note: \code{E} is reserved for residual (or environmental) variance and can be missed in \code{v}; it is considered to be an identify matrix if it is missing.

\code{v} can be provided as a single matrix, representing a variance component other than \code{E}.
}
  \item{initpar}{
Optional initial parameter values. When provided, \code{\link{optim}} will be called for optimization, which may take time but is good for checking of the result (see details for more).
}
  \item{nit}{
Maximum number of iterations for optimization. Ignored if there are not more than two variance components.
}
  \item{method}{
Either maximum likelihood (ML) or restricted maximum likelihood (REML).
}
  \item{control}{
A list of control parameters to be passed to \code{\link{optim}}.
}
  \item{hessian}{
Logical. Should a numerically differentiated Hessian matrix be returned?
}
}
\details{
The optimization function \code{\link{optim}} is adopted in the above function to estimate the parameters and maximum likelihood. Several optimization methods are available for the optimization algorithm in \code{\link{optim}}, but we recommend "Nelder-Mead" for the sake of stability. Alternatively, one may choose other options, e.g., "BFGS" to initialize and speed up the estimation procedure and then the procedure will automatically turn to "Nelder-Mead" for final results. If there is only one variance component (other than \code{E}), \code{\link{optimize}} will be used for optimization unless \code{initpar} is provided.

Normality is assumed for the random effects. Input data should be free of missing values.
}
\note{
   Hessian matrix, if requested, pertains to -log-likelihood function.
}
\value{
\item{par}{estimates of the model parameters.}
\item{value}{log-likelihood of the model.}
\item{y}{y used.}
\item{x}{associated with x used.}
\item{v}{variance component matrices v used.}
\item{\dots}{other information.}
}

\seealso{
\code{\link{optim}} and \code{\link{rem}}.
}
\examples{
data(miscEx)

\dontrun{
# no sex effect
pheno<- pdatF8[!is.na(pdatF8$bwt) & !is.na(pdatF8$sex),]
ii<- match(rownames(pheno), rownames(gmF8$AA))
v<- list(A=gmF8$AA[ii,ii], D=gmF8$DD[ii,ii])

o<- estVC(y=pheno$bwt, v=v)
o

# sex as fixed effect
fo<- estVC(y=pheno$bwt, x=pheno$sex, v=v)
fo
2*(fo$value-o$value) # log-likelihood test statistic

# sex as random effect
SM<- rem(~sex, data=pheno)
ro<- estVC(y=pheno$bwt, v=c(v,list(Sex=SM$sex)))
ro
2*(ro$value-o$value) # log-likelihood test statistic
}
}


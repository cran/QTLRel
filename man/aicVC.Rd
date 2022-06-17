\name{aicVC}
\alias{aicVC}
\title{
AIC Model Selection
}
\description{
Select genetic variance components via Akaike's information criterion (AIC).
}
\usage{
aicVC(y, x, v = list(E=diag(length(y))), initpar, k = 2, init = 1, keep = 1,
   direction = c("forward", "backward"), nit = 25, msg = FALSE,
   control = list(), hessian = FALSE)
}
\arguments{
  \item{y}{
A numeric vector or a numeric matrix of one column (representing a phenotype for instance).
}
  \item{x}{
A data frame or matrix, representing covariates if not missing.
}
  \item{v}{
A list of variance components of interest. Note: \code{E} is reserved for residual (or environmental) variance and can be missed in \code{v}; it is considered to be an identify matrix if it is specified.

\code{v} can be provided as a single matrix.
}
  \item{initpar}{
Optional initial parameter values.
}
  \item{k}{
Penalty on a parameter. The selection criterion is the known "AIC" if \code{k = 2} and is "BIC" if \code{k = log(n)} where "n" is the sample size.
}
  \item{init}{
Indicates which variance components for the initial model. By default, \code{E} is included if it is missing in \code{v}.
}
 \item{keep}{
Indicator of which variance components should be forced into the final model. By default, \code{E} is kept in the final model if it is not specified in \code{v}.
}
  \item{direction}{
The mode of search. Either "forward" or "backward" with default "forward".
}
  \item{nit}{
Maximum number of iterations for optimization. Ignored if there are not more than two variance components.
}
  \item{msg}{
A logical variable. True if one wants to track the process for monitoring purpose.
}
  \item{control}{
A list of control parameters to be passed to \code{\link{optim}}.
}
  \item{hessian}{
Logical. Should a numerically differentiated Hessian matrix be returned?
}
}
\details{
In genome-wide association studies (GWAS), random effects are usually added to a model to account for polygenic variation. Abney et al (2000) showed that five variance components including the most interesting additive and dominance variance components are potentially induced by polygenes. The above function is intended for selecting variance components that contribute "most" to a quantitative trait.

Function \code{\link{estVC}} is called by the above function to estimate the parameters and maximum likelihood in each model. Refer to \code{\link{estVC}} for more information.
}
\value{
\item{aic}{AIC of the final model.}
\item{model}{Gives parameter estimates, log-likihood, and other information.}
\item{lik}{Log-likelihood of the model selected at each intermediate step.}
\item{trace}{Indicates which variance components were selected at each intermediate step.}
}
\seealso{
\code{\link{estVC}} for more information.
}
\examples{
data(miscEx)

\dontrun{
# forward selection
# any variance component will be selected
# if AIC improve by 1e-5 or larger
pheno<- pdatF8[!is.na(pdatF8$bwt) & !is.na(pdatF8$sex),]
ii<- match(rownames(pheno), rownames(gmF8$AA))
v<- list(A=gmF8$AA[ii,ii], D=gmF8$DD[ii,ii])

o<- aicVC(y=pheno$bwt, x=pheno$sex, k=0, v=v, msg=TRUE)
o

# forward selection
of<- aicVC(y=pheno$bwt, x=pheno$sex, v=v, k=1/2,
	direction="for", msg=TRUE)
of

# backward elimination
ob<- aicVC(y=pheno$bwt, x=pheno$sex, v=v, k=1/2, init=1:2,
	direction="back", msg=TRUE)
ob
}
}


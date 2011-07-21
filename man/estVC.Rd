\name{estVC}
\alias{estVC}
\title{
Estimate Variance Component Parameters
}
\description{
Estimate model parameters for covariates, genetic variance components and residual effect.
}
\usage{
estVC(y,x,v=vector("list",6),initpar,nit=25,
   method=c("Nelder-Mead","BFGS","CG","SANN"),
   control=list(),hessian=FALSE)
}
\arguments{
  \item{y}{
a numeric vector or a numeric matrix of one column (representing a phenotype for instance).
}
  \item{x}{
a data frame or matrix, representing covariates if not missing.
}
  \item{v}{
a list of variance components (AA, DD, HH, AD, MH, EE,...), where "AA" and "DD" are respectively additive and dominance genetic matrices, "HH", "AD" and "MH" are other genetic matrices that one may be interested in (see \code{\link{aicVC}}), "EE" is the residual matrix that is usually assumed to be an identity matrix, and "..." are other random components of interest. If a genetic component is not considered, it should be set to NULL.
}
  \item{initpar}{
optional initial parameter values.
}
  \item{nit}{
number of iterations to call \code{\link{optim}} for optimization.
}
  \item{method}{
the optimization method to be used. See \code{\link{optim}} for details.
}
  \item{control}{
a list of control parameters to be passed to \code{\link{optim}}.
}
  \item{hessian}{
logical. Should a numerically differentiated Hessian matrix be returned?
}
}
\details{
The optimization function \code{\link{optim}} is adopted in the above function to estimate the parameters and maximum likelihood. Several optimization methods are available for the optimization algorithm in \code{\link{optim}}, but we recommend "Nelder-Mead" for the sake of stability. Alternatively, one may choose other options, e.g., "BFGS" to initialize and speed up the estimation procedure and then the procedure will automatically turn to "Nelder-Mead" for final results.

Normality is assumed for the random effects. Input data should be free of missing values.
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
x<- sample(c("M","F"),100,replace=TRUE)
y<- rnorm(100) + (x=="M")*1.0
v<- cov(matrix(rnorm(10^6),ncol=100))

\dontrun{
# no sex effect
o<- estVC(y, v = list(AA=v,DD=NULL,HH=NULL,AD=NULL,
   MH=NULL,EE=diag(100)),method="BFGS")
o

# sex as fixed effect
fo<- estVC(y, x, v = list(AA=v,DD=NULL,HH=NULL,AD=NULL,
   MH=NULL,EE=diag(100)),method="BFGS")
fo
2*(fo$value-o$value) # log-likelihood test statistic

# sex as random effect
SM<- rem(~x,data.frame(x=x))
ro<- estVC(y, v = list(AA=v,DD=NULL,HH=NULL,AD=NULL,
   MH=NULL,SE=SM$x,EE=diag(100)),method="BFGS")
ro
2*(ro$value-o$value) # log-likelihood test statistic
}
}


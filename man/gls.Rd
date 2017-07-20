\name{gls}
\alias{gls}
\title{
Generalized Least Squares Estimates
}
\description{
Obtain estimates using generalized least squares (gls).
}
\usage{
gls(formula, data, vc = NULL, test=c("none","F"))
}
\arguments{
  \item{formula}{
An object of class "formula": a symbolic description of the model to be fitted.
}
  \item{data}{
An data frame containing the variables in the model.
}
  \item{vc}{
An object from \code{\link{estVC}} or \code{\link{aicVC}} or an estimated variance-covariance matrix induced by relatedness and environment if not NULL.
}
  \item{test}{
Wheter F-test is performed.
}
}
\value{
A matrix with columns: "Estimate", "Std. Error", "t value" and "Pr(>|t|)", or an ANOVA table if F-test is requested.
}

\seealso{
\code{\link{lm}}.
}


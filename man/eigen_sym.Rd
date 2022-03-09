\name{eigen.sym}
\alias{eigen.sym}
\title{
Spectral decomposition of a matrix
}
\description{
Computes eigenvalues and eigenvectors of real symmetric matrices.
}
\usage{
eigen.sym(x)
}
\arguments{
  \item{x}{
A real symmetric matrix.
}
}
\details{
This is to use the LAPACK routine 'DSYEVR' to perform spectral decomposition.
}
\value{
\item{values}{a vector containing the eigenvalues of \code{x}, sorted in decreasing order.}
\item{vectors}{a matrix whose columns contain the eigenvectors of \code{x}, corresponding to eigenvalues.}
}
\note{
Warning: symmetry is not checked by the program!
}
\seealso{
\code{\link{eigen}} for more information.
}


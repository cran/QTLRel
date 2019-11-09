\name{pedRecode}
\alias{pedRecode}
\title{
Recode a Pedigree
}
\description{
Prepare a pedigree in a format that is suitable for certain functions
}
\usage{
pedRecode(ped, ids, all = TRUE, msg = TRUE)
}
\arguments{
  \item{ped}{
A pedigree, which is a data frame (id, father/sire, mother/dam, ...). If "sex" is a component, male should be "M", "Male" or 1, and female should be "F", "Female" or 2 (other than 0 and 1). If given, "generation" can be numeric 0, 1, 2, ... or non-numeric "F0", "F1", "F2", ..., which should be in an increasing order. Note: 0 is reserved for unknown \code{father}, \code{mother} or \code{sex}. If a father/mother is an inbred founder, its ID should be tagged by character 'i' (e.g. 1i, 2i, etc.).
}
  \item{ids}{
If given, only individuals with \code{ids} and their ancestors are kept in the recoded pedigree.
}
  \item{all}{
If false, fathers and mothers with no parents are treated as unknown.
}
  \item{msg}{
If false, messages are suppressed.
}
}
\details{
This function is used in \code{cic}, and it can be used for error checking with respect to sex and generation if sex and/or generation information is available. The actual values of generation can be anything but should correspond to the true order of generation; otherwise, \code{cic} may fail or we may get incorrect results. Information except \code{id}, \code{father} and \code{mother} is optional.
}
\value{
A recoded pedigree.
}
\seealso{
\code{\link{cic}}.
}

\examples{
data(miscEx)

pedF8[1:10,]
pedR<- pedRecode(pedF8)
pedR[1:10,]
dim(pedR)
pedR<- pedRecode(pedF8, ids=pedF8$id[pedF8$gener=="F8"])
dim(pedR)
}


\name{pedRecode}
\alias{pedRecode}
\title{
Recode a Pedigree
}
\description{
Prepare a pedigree in a format that is suitable for certain functions
}
\usage{
pedRecode(ped,ids)
}
\arguments{
  \item{ped}{
a pedigree, which is a data frame (id, sire, dam, ...). If "sex" is a component, male should be "M", "Male" or 1, and female should be "F", "Female" or 2 (other than 0 and 1). If given, "generation" can be numeric 0, 1, 2, ... or non-numeric "F0", "F1", "F2", etc. 0 is reserved for unknown \code{sire} or \code{dam}.
}
  \item{ids}{
If given, only individuals with \code{ids} and their ancestors are kept in the recoded pedigree.
}
}
\value{
A recoded pedigree.
}
\examples{
data(miscEx)

ped[1:10,]
pedR<- pedRecode(ped)
pedR[1:10,]
dim(pedR)
pedR<- pedRecode(ped,ids=ped$id[ped$gener=="F34"])
dim(pedR)
}


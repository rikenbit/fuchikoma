\docType{methods}
\name{CatKernel}
\alias{CatKernel}
\title{
Categorical Kernel function
}
\description{
CatKernel calculate gram matrix by using label information (e.g., cell label or sample label such as c(1,1,1,2,2,2,3,3)).
}
\usage{
CatKernel(label, type = c("simple", "two", "one_vs_rest", "each"))
}
\arguments{
  \item{label}{
Label vector descrbing which cells are sample cell type and which cells are different (e.g., c(1,1,1,2,2,2,3,3) (default:simple)).
}
  \item{type}{
Categorical kernel functions. simple could be specified even the number of class is two or more. When the number of class is two, only two parameter is accessible. When the number of class is three or more, one_vs_rest and each is accessible.
}
}
\value{
A data frame is returned.
}
\references{
Le Song et al. (2007) Gene selection via the BAHSIC family of algorithms, Bioinformatics, 23(13), i490-i498

Y-h Taguchi et al. (2015) Principal component analysis-based unsupervised feature extraction applied to in silico drug discovery for posttraumatic stress disorder-mediated heart disease, BMC Bioinformatics, 16(139)
}
\author{
\packageAuthor{fuchikoma}

Maintainer: \packageMaintainer{fuchikoma}
}
\examples{
data(label.MARS)
L <- CatKernel(label.MARS)
dim(L)
}
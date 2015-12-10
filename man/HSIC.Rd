\docType{methods}
\name{HSIC}
\alias{HSIC}
\title{
Hilbert-Schmidt Independence Criteria (HSIC)
}
\description{
Inner function of FUCHIKOMA. HSIC calculated aginst two gram-matrix means the independence of two kernel spaces. The higher HSIC value is, the more two space are dependent. When p.value set TRUE, p-value of HSIC is also calculated by moment matching to a gamma distribution.
}
\usage{
HSIC(K, L, shrink=FALSE, type=c("gamma", "permutation"), n.perm=100)
}
\arguments{
  \item{K}{First gram matrix}
  \item{L}{Second gram matrix}
  \item{shrink}{Shrinkage option when the sample size is insufficient}
  \item{type}{When type is gamma, p-value is calculated under the assumption of gamma distribution. When type is permutation, p-value is calculated by random permutation test}
  \item{n.perm}{Number of permutation (default:100)}
 }
\value{
\item{HSIC}{The value of HSIC}
\item{Pval}{P-value of HSIC. The value is accessible only when p.value is specified as TRUE}
}
\references{
Le Song et al. (2007) Gene selection via the BAHSIC family of algorithms, Bioinformatics, 23(13), i490-i498

Y-h Taguchi et al. (2015) Principal component analysis-based unsupervised feature extraction applied to in silico drug discovery for posttraumatic stress disorder-mediated heart disease, BMC Bioinformatics, 16(139)

Arthur Gretton et al. (2007) A Kernel Statistical Test of Independence, NIPS 21

Aaditya Ramdas et al. (2015) Nonparametric Independence Testing for Small Sample Sizes, IJCAI-15
}
\author{
\packageAuthor{FUCHIKOMA}

Maintainer: \packageMaintainer{FUCHIKOMA}
}

\examples{
K <- matrix(runif(100), nrow=10)
L <- matrix(runif(100), nrow=10)
HSIC(K, L)
}
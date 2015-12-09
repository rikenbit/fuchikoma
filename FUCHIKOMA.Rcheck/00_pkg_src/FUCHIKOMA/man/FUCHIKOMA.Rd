\docType{methods}
\name{FUCHIKOMA}
\alias{FUCHIKOMA}
\title{
\packageTitle{FUCHIKOMA}
}
\description{
\packageDESCRIPTION{FUCHIKOMA}
}
\usage{
FUCHIKOMA(data, mode = c("Supervised", "Unsupervised"), Comp = FALSE, label = FALSE, cat.type = FALSE, n.eigs = 10, algorithm = c("brute", "song"), per.rej = 10, threshold = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{A data matrix, in which the row means genes and column means cells or sample.}
\item{mode}{When Supervised is specified, FUCHIKOMA uses label paramter. When Unsupervised is specified, FUCHIKOMA uses Comp parameter for specifying which diffusion components should be used.}
\item{Comp}{When mode is specified as Unsupervised, Comp must be specified such as c(1,2).}
\item{label}{When mode is specified as Supervised, label must be specified such as c(1,1,1,2,2,2,3,3)}
\item{cat.type}{Type of categorical kernel of \code{CatKernel}}
\item{n.eigs}{Number of eigenvectors/values to return (default: 20)}
\item{algorithm}{brute means single gene rejection strategies and song means fixed percent of genes rejection strategies in each iteration step of FUCHIKOMA.
}
\item{per.rej}{When algorithm is specified as song, per.rej must be specified such as 20 (default:10)}
\item{threshold}{In each iteration step, if the difference of HSIC in the step and max value of previous HSICs is lower than threshold, iteration will be halted (default: 0.01).
}
}
\value{
\item{DEGs.HSICs}{Differentially expressed genes and HSIC values}
\item{DEGs.Pvals}{Pvalues of HSICs}
\item{All.HSICs}{HSICs value of all genes}
\item{Rej.order}{The order of rejecting of genes in each iteration step}
}
\references{
Laleh Haghverdi et al. (2015) Diffusion maps for high-dimensional single-cell analysis of differentiation data. Bioinformatics, 31(18), 2989-2998

Le Song et al. (2007) Gene selection via the BAHSIC family of algorithms, Bioinformatics, 23(13), i490-i498

Y-h Taguchi et al. (2015) Principal component analysis-based unsupervised feature extraction applied to in silico drug discovery for posttraumatic stress disorder-mediated heart disease, BMC Bioinformatics, 16(139)
}
\author{
\packageAuthor{FUCHIKOMA}

Maintainer: \packageMaintainer{FUCHIKOMA}
}

\seealso{
\code{\link[destiny]{DiffusionMap}}
}
\examples{
data(MARS)
res <- FUCHIKOMA(data=MARS, mode="Supervised", label=label.MARS, type="one_vs_rest", n.eigs=10, algorithm="song", per.rej=10, threshold=0.01)
}
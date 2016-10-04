\docType{methods}
\name{fuchikoma}
\alias{fuchikoma}
\title{
\packageTitle{fuchikoma}
}
\description{
\packageDESCRIPTION{fuchikoma}
}
\usage{
fuchikoma(data, cores=NULL, mode = c("Supervised", "Unsupervised", "Mix", "tSNE"), weight=c(0.5,0.5), Comp = NULL, label = FALSE, cat.type = c("simple", "one_vs_rest", "each", "two"), n.eigs = 10, algorithm = c("song", "brute"), per.rej = 10, threshold = 0.01, verbose=FALSE, dropout=10, sigma=15, perplexity = 30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{A data matrix, in which the row means genes and column means cells or sample.}
\item{cores}{The number of using mathine cores. (default: automaticaly detected number by doParallel package)}
\item{mode}{When Supervised is specified, fuchikoma uses label paramter. When Unsupervised is specified, fuchikoma uses Comp parameter for specifying which diffusion components should be used. When Mixed is specified, fuchikoma uses both mode (default:Supervised)}
\item{weight}{Weight for integrating two gram matrix, when mode is specified as Mix (default:0.5 vs 0.5)}
\item{Comp}{When mode is specified as Unsupervised, Comp must be specified such as c(1,2). (default:FALSE)}
\item{label}{When mode is specified as Supervised, label must be specified such as c(1,1,1,2,2,2,3,3) (default:FALSE)}
\item{cat.type}{Type of categorical kernel of \code{CatKernel} (default: simple)}
\item{n.eigs}{Number of eigenvectors/values to return (default: 20)}
\item{algorithm}{brute means single gene rejection strategies and song means fixed percent of genes rejection strategies in each iteration step of fuchikoma.
 (default:song)}
\item{per.rej}{When algorithm is specified as song, per.rej must be specified such as 20 (default:10)}
\item{threshold}{In each iteration step, if the difference of HSIC in the step and max value of previous HSICs is lower than threshold, iteration will be halted (default: 0.01).
}
\item{verbose}{verbose option (default: FALSE).}
\item{dropout}{Threshold when the remaining gene is few (default: 10).}
\item{sigma}{Parameter used in DiffusionMap of destiny (default: 15).}
\item{perplexity}{Parameter of tSNE (default: 30).}
}
\value{
\item{DEGs.HSICs}{Differentially expressed genes (DEGs) and HSIC values}
\item{DEGs.Pvals}{Pvalues of DEGs}
\item{All.HSICs}{HSICs value of all genes}
\item{All.Pvals}{Pvalues of all genes}
}
\references{
Laleh Haghverdi et al. (2015) Diffusion maps for high-dimensional single-cell analysis of differentiation data. Bioinformatics, 31(18), 2989-2998

Le Song et al. (2007) Gene selection via the BAHSIC family of algorithms, Bioinformatics, 23(13), i490-i498

Y-h Taguchi et al. (2015) Principal component analysis-based unsupervised feature extraction applied to in silico drug discovery for posttraumatic stress disorder-mediated heart disease, BMC Bioinformatics, 16(139)

Aaditya Ramdas et al. (2015) Nonparametric Independence Testing for Small Sample Sizes, IJCAI-15
}
\author{
\packageAuthor{fuchikoma}

Maintainer: \packageMaintainer{fuchikoma}
}

\seealso{
\code{\link[destiny]{DiffusionMap}}
}
\examples{
# three cell data
CellA <- data.frame(matrix(rnorm(50*20), nrow=50, ncol=20))
CellB <- data.frame(matrix(rnorm(50*20), nrow=50, ncol=20))

# DEGs definition
CellA[1:10, ] <- CellA[1:10, ] + 10 * matrix(runif(10*20), nrow=10, ncol=20)
CellB[11:20, ] <- CellB[1:10, ] + 10 * matrix(runif(10*20), nrow=10, ncol=20)
# testdata
testdata2 <- data.frame(CellA, CellB)
colnames(testdata2) <- c(paste0("CellA_", 1:20), paste0("CellB_", 1:20))
rownames(testdata2) <- paste0("Gene", 1:nrow(testdata2))

# label
label2 <- c(rep(1, 20), rep(2, 20))

res <- fuchikoma(data=testdata2, cores=1, label=label2)
}
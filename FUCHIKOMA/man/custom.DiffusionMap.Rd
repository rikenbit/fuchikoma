\docType{methods}
\name{custom.DiffusionMap}
\alias{custom.DiffusionMap}
\title{
Customized Diffusion Map of destiny for extracting Gram Matrix
}
\description{
Inner function of FUCHIKOMA, which returns the gram matrix of diffusion map.
}
\usage{
custom.DiffusionMap(data, sigma = NULL, k = find.dm.k(nrow(data) - 1L), n.eigs = min(20L, nrow(data) - 2L), density.norm = TRUE, ..., distance = c("euclidean", "cosine"), censor.val = NULL, censor.range = NULL, missing.range = NULL, vars = NULL, verbose = !is.null(censor.range), .debug.env = NULL)
}
\arguments{
\item{data}{Expression data to be analyzed. Provide \code{vars} to select specific columns other than the default: all double value columns}

\item{sigma}{Diffusion scale parameter of the Gaussian kernel. Either a number or a \link{Sigmas} object.
(Optional. default: will be calculated using \link{find.sigmas})
A larger sigma might be necessary if the eigenvalues can not be found because of a singularity in the matrix}

\item{k}{Number of nearest neighbors to consider (default: a guess betweeen 100 and \eqn{n - 1})
\code{NULL} or \code{NA} are also interpreted as \code{n - 1L}.}

\item{n.eigs}{Number of eigenvectors/values to return (default: 20)}

\item{density.norm}{logical. If TRUE, use density normalisation}

\item{...}{All parameter after this are optional and have to be specified by name}

\item{distance}{Distance measurement method. Euclidean distance (default) or cosine distance (\eqn{1-corr(c_1, c_2)}).}

\item{censor.val}{Value regarded as uncertain. Either a single value or one for every dimension (Optional, default: CENSOR.VAL)}

\item{censor.range}{Uncertainity range for censoring (Optional, default: none). A length-2-vector of certainty range start and end. TODO: also allow \eqn{2\times G} matrix}

\item{missing.range}{Whole data range for missing value model. Has to be specified if NAs are in the data}

\item{vars}{Variables (columns) of the data to use. Specifying NULL will select all columns (default: All floating point value columns)}

\item{verbose}{Show a progressbar and other progress information (default: do it if censoring is enabled)}

\item{.debug.env}{If supplied, the rotated transition probaility matrix \code{M} and the scaling rotation \code{D.rot} is stored in it prior to eigen decomposition.}
}
\value{
A DiffusionMap object:
}
\references{
Laleh Haghverdi et al. (2015) Diffusion maps for high-dimensional single-cell analysis of differentiation data. Bioinformatics, 31(18), 2989-2998
}
\author{
\packageAuthor{FUCHIKOMA}

Maintainer: \packageMaintainer{FUCHIKOMA}
}

\seealso{
\code{\link[destiny]{DiffusionMap}}
}
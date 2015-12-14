\name{result.dmap.MARS}
\alias{result.dmap.MARS}
\docType{data}
\title{
The result of Diffusion Map against MARS-Seq data
}
\description{
The result of dmap against MARS.

The data has four different cell types; B cells (B1-B48), Dendritic cells (DC1-DC89), Monocytes (Mono1-Mono46), and Natural Killer cells (NK1-NK45).
}
\usage{data(result.dmap.MARS)}
\source{
\url{https://www.sciencemag.org/content/343/6172/776?related-urls=yes&legid=sci;343/6172/776}

\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54006}
}
\author{
\packageAuthor{fuchikoma}

Maintainer: \packageMaintainer{fuchikoma}
}
\references{
Diego Adhemar Jaitin, Ephraim Kenigsberg, Hadas Keren-Shaul, Naama Elefant, Franziska Paul, Irina Zaretsky, Alexander Mildner, Nadav Cohen, Steffen Jung, Amos Tanay, Ido Amit (2014) Massively Parallel Single-Cell RNA-Seq for Marker-Free Decomposition of Tissues into Cell Types. \emph{Science}, \bold{343 (6172)}: 776-779
}
\examples{
data(MARS)
data(label.MARS)
data(result.pca.MARS)
data(result.tsne.MARS)
data(result.dmap.MARS)
pairs(result.pca.MARS, col=label.MARS, main="MARS-Seq (PCA)", pch=16)
plot(result.tsne.MARS, col=label.MARS, main="MARS-Seq (t-SNE)", pch=16)
pairs(result.dmap.MARS, col=label.MARS, main="MARS-Seq (Diffusion Map)", pch=16)
}
\keyword{datasets}
\name{result.tsne.MARS}
\alias{result.tsne.MARS}
\docType{data}
\title{
The result of tsne against MARS-Seq data
}
\description{
The result of t-SNE against MARS.

The data has four different cell types; B cells (B1-B48), Dendritic cells (DC1-DC89), Monocytes (Mono1-Mono46), and Natural Killer cells (NK1-NK45).
}
\usage{data(result.tsne.MARS)}
\source{
\url{https://www.sciencemag.org/content/343/6172/776?related-urls=yes&legid=sci;343/6172/776}

\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54006}
}
\author{
\packageAuthor{fuchikoma}

Maintainer: \packageMaintainer{fuchikoma}
}
\references{
L. J. P. van der Maaten et al. (2008) Visualizing High-Dimensional Data Using t-SNE. Journal of Machine Learning Research, 9(Nov), 2579-2605
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
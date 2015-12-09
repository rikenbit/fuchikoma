pkgname <- "FUCHIKOMA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('FUCHIKOMA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("CatKernel")
### * CatKernel

flush(stderr()); flush(stdout())

### Name: CatKernel
### Title: Categorical Kernel function
### Aliases: CatKernel

### ** Examples

data(MARS)
L <- CatKernel(label.MARS, type="one_vs_rest")
dim(L)



cleanEx()
nameEx("FUCHIKOMA")
### * FUCHIKOMA

flush(stderr()); flush(stdout())

### Name: FUCHIKOMA
### Title: Detection of Differentially Expressed Genes in one of multiple
###   clusters using BAHSIC and Diffusion Map
### Aliases: FUCHIKOMA

### ** Examples

data(MARS)
res <- FUCHIKOMA(data=MARS, mode="Supervised", label=label.MARS, type="one_vs_rest", n.eigs=10, algorithm="song", per.rej=10, threshold=0.01)



cleanEx()
nameEx("HSIC")
### * HSIC

flush(stderr()); flush(stdout())

### Name: HSIC
### Title: Hilbert-Schmidt Independence Criteria (HSIC)
### Aliases: HSIC

### ** Examples

K <- matrix(runif(100), nrow=10)
L <- matrix(runif(100), nrow=10)
HSIC(K, L, p.value)



cleanEx()
nameEx("MARS")
### * MARS

flush(stderr()); flush(stdout())

### Name: MARS
### Title: Count data of MARS-Seq containing four different cell types.
### Aliases: MARS
### Keywords: datasets

### ** Examples

data(MARS)
dim(MARS)
pairs(result.pca.MARS$rotation[,1:10], col=label.MARS, main="MARS-Seq (PCA)")
pairs(result.destiny.MARS@eigenvectors[,1:10], col=label.MARS, main="MARS-Seq (Diffusion Map)")



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

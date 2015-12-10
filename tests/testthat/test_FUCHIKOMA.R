######### testdata ########
# three cell data
CellA <- data.frame(matrix(rnorm(50*20), nrow=50, ncol=20))
CellB <- data.frame(matrix(rnorm(50*20), nrow=50, ncol=20))
CellC <- data.frame(matrix(rnorm(50*20), nrow=50, ncol=20))

# DEGs definition
CellA[1:10, ] <- CellA[1:10, ] + 10 * matrix(runif(10*20), nrow=10, ncol=20)
CellB[11:20, ] <- CellB[1:10, ] + 10 * matrix(runif(10*20), nrow=10, ncol=20)
CellC[21:30, ] <- CellC[1:10, ] + 10 * matrix(runif(10*20), nrow=10, ncol=20)

# testdata
testdata2 <- data.frame(CellA, CellB)
testdata3 <- data.frame(CellA, CellB, CellC)
colnames(testdata2) <- c(paste0("CellA_", 1:20), paste0("CellB_", 1:20))
colnames(testdata3) <- c(paste0("CellA_", 1:20), paste0("CellB_", 1:20), paste0("CellC_", 1:20))
rownames(testdata2) <- paste0("Gene", 1:nrow(testdata2))
rownames(testdata3) <- paste0("Gene", 1:nrow(testdata3))

# label
label <- c(rep(1, 20), rep(2, 20), rep(3, 20))

######### All combination ########
##### 2 group
### Supervised

### UnSupervised


##### 3 group
### Supervised

### UnSupervised


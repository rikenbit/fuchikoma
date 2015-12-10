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
testdata <- data.frame(CellA, CellB, CellC)
colnames(testdata) <- c(paste0("CellA_", 1:20), paste0("CellB_", 1:20), paste0("CellC_", 1:20))
rownames(testdata) <- paste0("Gene", 1:nrow(testdata))

# label
label <- c(rep(1, 20), rep(2, 20), rep(3, 20))

######### Normal run ########
res1 <- FUCHIKOMA(data=testdata, mode="Supervised", label=label, cat.type="one_vs_rest", n.eigs=10, algorithm="song", verbose=TRUE)

expect_equivalent(2, 1)

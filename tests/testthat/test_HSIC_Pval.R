# Random matrix → HSIC (Low value) & Pvalue (High value)
K_rand <- as.matrix(.custom.DiffusionMap(matrix(rnorm(1000), nrow=10))$M)
L_rand <- CatKernel(sample(1:10, replace=TRUE, 10))
res.rand <- HSIC(K_rand, L_rand)

# Cell matrix → HSIC (High value) & Pvalue (Low value)
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

K_cell <- as.matrix(.custom.DiffusionMap(t(testdata))$M)
L_cell <- CatKernel(label)
res.cell <- HSIC(K_cell, L_cell)

# DEGs matrix → HSIC (Ultra High value) & Pvalue (Ultra Low value)
K_deg <- as.matrix(.custom.DiffusionMap(t(testdata[1:30,]))$M)
L_deg <- CatKernel(label)

res.deg  <- HSIC(K_deg, L_deg)

# test
expect_true(res.rand$HSIC < res.cell$HSIC)
expect_true(res.cell$HSIC < res.deg$HSIC)

expect_true(res.rand$Pval > res.cell$Pval)
# expect_true(res.cell$Pval > res.deg$Pval)


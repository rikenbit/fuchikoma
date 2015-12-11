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
label2 <- c(rep(1, 20), rep(2, 20))
label3 <- c(rep(1, 20), rep(2, 20), rep(3, 20))

######### All combination ########
############# 2 group ############
### Supervised
res2_1 <- FUCHIKOMA(data=testdata2, label=label2)
res2_2 <- FUCHIKOMA(data=testdata2, label=label2, cat.type="two")

expect_warning(expect_error(res2_3 <- FUCHIKOMA(data=testdata2, label=label2, cat.type="one_vs_rest")))
expect_warning(expect_error(res2_4 <- FUCHIKOMA(data=testdata2, label=label2, cat.type="each")))

### UnSupervised
res2_5 <- FUCHIKOMA(data=testdata2, mode="Unsupervised", Comp=c(1,2,3))
expect_error(res2_6 <- FUCHIKOMA(data=testdata2, mode="Unsupervised", Comp=1:100))

############# 3 group ############
### Supervised
res3_1 <- FUCHIKOMA(data=testdata3, label=label3)
expect_warning(expect_error(res3_2 <- FUCHIKOMA(data=testdata3, label=label3, cat.type="two")))

res3_3 <- FUCHIKOMA(data=testdata3, label=label3, cat.type="one_vs_rest")
res3_4 <- FUCHIKOMA(data=testdata3, label=label3, cat.type="each")

### UnSupervised
res3_5 <- FUCHIKOMA(data=testdata3, mode="Unsupervised", Comp=c(1,2,3))
expect_error(res3_6 <- FUCHIKOMA(data=testdata3, mode="Unsupervised", Comp=1:100))

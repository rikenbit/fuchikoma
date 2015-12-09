# # パッケージ読み込み
# source("http://bioconductor.org/biocLite.R")
# biocLite("destiny")
# install.packages("foreach")
# install.packages("doParallel")

# library("destiny")
# library("foreach")
# library("doParallel")

####################### 試しに実行 #####################
# 標準正規分布
CellA <- data.frame(matrix(rnorm(100*20), nrow=100, ncol=20))
CellB <- data.frame(matrix(rnorm(100*20), nrow=100, ncol=20))
CellC <- data.frame(matrix(rnorm(100*20), nrow=100, ncol=20))

# DEGを指定（もっと最適な作り方を調べる必要あり）
CellA[1:10, ] <- CellA[1:10, ] + 10 * matrix(runif(10*20), nrow=10, ncol=20)
CellB[11:20, ] <- CellB[1:10, ] + 10 * matrix(runif(10*20), nrow=10, ncol=20)
CellC[21:30, ] <- CellC[1:10, ] + 10 * matrix(runif(10*20), nrow=10, ncol=20)

testdata <- data.frame(CellA, CellB, CellC)

colnames(testdata) <- c(paste0("CellA_", 1:20), paste0("CellB_", 1:20), paste0("CellC_", 1:20))
rownames(testdata) <- paste0("Gene", 1:nrow(testdata))

# ラベル
label <- c(rep(1, 20), rep(2, 20), rep(3, 20))

# # 普通のPCA
# res.pca <- prcomp(testdata)
# pairs(res.pca$rotation[, 1:20], col=label)



# ######## Diffusion Mapを実行（第10主成分までを見る） ####
# # オブジェクト化
# testdata.obj <- as.ExpressionSet(as.data.frame(t(testdata)))

# dif1 <- DiffusionMap(testdata.obj, n.eigs=10)
# pairs(eigenvectors(dif1), col=label)

# # 改造Ver Diffusion Map
# dif2 <- custom.DiffusionMap(testdata.obj, n.eigs=10)
# pairs(dif2$eigenvectors, col=label)



################### FUCHIKOMA実行 ###################
# 教師あり（クラスラベルを与える）
result1 <- FUCHIKOMA(data=testdata, mode="Supervised", label=label, cat.type="one_vs_rest", n.eigs=10, algorithm="song", per.rej=10, threshold=0.01)

result1
plot(result1$All.HSICs)

# result2 <- FUCHIKOMA(data=testdata, mode="Supervised", label=label, cat.type="each", n.eigs=10, algorithm="song", per.rej=5, threshold=0.01)

# result2
# plot(result2$All.HSICs)

# result3 <- FUCHIKOMA(data=testdata, mode="Supervised", label=label, cat.type="simple", n.eigs=10, algorithm="song", per.rej=5, threshold=0.01)

# result3
# plot(result3$All.HSICs)


# # 教師なし（指定した主成分を使う）
# result4 <- FUCHIKOMA(data=testdata, mode="Unsupervised", Comp=c(1,2,3), n.eigs=10, algorithm="song", per.rej=5, threshold=0.01)

# result4
# plot(result4$All.HSICs)


# ########## HSICのテスト1 ############
# ### ランダムな行列で、HSICが0に近く、###
# #### Pvalueが1に近くなるかどうか ######
# ####################################
# HSICs <- c()
# Pvals <- c()
# sameHSICs <- c()
# samePvals <- c()
# for(i in 1:100){
#     print(i)
#     A <- as.matrix(custom.DiffusionMap(matrix(rnorm(1000), nrow=10))$M)
#     B <- CatKernel(sample(1:10, replace=TRUE, 10), cat.type="one_vs_rest")

#     HSICs[i] <- HSIC(A, B)$HSIC
#     sameHSICs[i] <- HSIC(A, B)$HSIC
#     Pvals[i] <- HSIC(A, B, p.value=TRUE)$Pval
#     samePvals[i] <- HSIC(A, B, p.value=TRUE)$Pval
# }
# HSICs
# sameHSICs
# Pvals
# samePvals

# ################ HSICのテスト2 ##################
# ### 細胞を模したデータで、HSICがある程度大きくなり ###
# ######## Pvalueが0により近くなるかどうか ##########
# ###############################################
# A <- as.matrix(custom.DiffusionMap(t(testdata))$M)
# B <- CatKernel(label, cat.type="one_vs_rest")
# HSIC(A, B)$HSIC
# HSIC(A, A)$HSIC
# HSIC(B, B)$HSIC
# HSIC(A, B, p.value=TRUE)$Pval
# HSIC(A, A, p.value=TRUE)$Pval
# HSIC(B, B, p.value=TRUE)$Pval


# ########### HSIC・FUCHIKOMAのテスト1 ############
# ############ FUCHIKOMAの各ステップで、 ###########
# ########### HSICの値が徐々に大きくなり、 ##########
# ############ Pvalueが徐々に0に近づき #############
# ############ DEGs数近辺で最大値をとるか ###########
# ###############################################




# ########### HSIC・FUCHIKOMAのテスト2 ############
# ######### DEGsだけを使って、FUCHI
# #
# ###############################################


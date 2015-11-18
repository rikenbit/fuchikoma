# パッケージ読み込み
library(kernlab)
library(destiny)

###################### 各種関数定義 #####################

# グラム行列を出力するように改造したDiffusionMap
custom.DiffusionMap <- function(data, sigma = NULL, k = find.dm.k(nrow(data) - 1L),
    n.eigs = min(20L, nrow(data) - 2L), density.norm = TRUE,
    ..., distance = c("euclidean", "cosine"), censor.val = NULL,
    censor.range = NULL, missing.range = NULL, vars = NULL, verbose = !is.null(censor.range),
    .debug.env = NULL)
{
    distance <- match.arg(distance)
    data.env <- new.env(parent = .GlobalEnv)
    data.env$data <- data
    # パッケージ名を指定
    data <- destiny:::extract.doublematrix(data, vars)
    imputed.data <- data
    if (any(is.na(imputed.data)))
        imputed.data <- as.matrix(hotdeck(data, imp_var = FALSE))
    n <- nrow(imputed.data)
    if (n <= n.eigs + 1L)
        stop(sprintf("Eigen decomposition not possible if n ≤ n.eigs+1 (And %s ≤ %s)",
            n, n.eigs + 1L))
    if (is.null(k) || is.na(k))
        k <- n - 1L
    if (k >= nrow(imputed.data))
        stop(sprintf("k has to be < nrow(data) (And %s ≥ nrow(data))",
            k))
    # パッケージ名を指定
    censor <- destiny:::test.censoring(censor.val, censor.range, data,
        missing.range)
    if (identical(distance, "cosine") && censor)
        stop("cosine distance only valid without censoring model")
    sigmas <- sigma
    if (is.null(sigmas)) {
        sigmas <- find.sigmas(imputed.data, censor.val = censor.val,
            censor.range = censor.range, missing.range = missing.range,
            vars = vars, verbose = verbose)
    }
    else if (!is(sigmas, "Sigmas")) {
        sigmas <- new("Sigmas", log.sigmas = NULL, dim.norms = NULL,
            optimal.sigma = sigma, optimal.idx = NULL, avrd.norms = NULL)
    }
    sigma <- optimal.sigma(sigmas)
    if (verbose) {
        cat("finding knns...")
        tic <- proc.time()
    }
    # パッケージ名を指定
    knn <- FNN:::get.knn(imputed.data, k, algorithm = "cover_tree") # <-
    if (verbose) {
        cat("...done. Time:\n")
        print(proc.time() - tic)
    }
    cb <- invisible
    if (verbose) {
        pb <- txtProgressBar(1, n, style = 3)
        cb <- function(i) setTxtProgressBar(pb, i)
        cat("Calculating transition probabilities...\n")
        tic <- proc.time()
    }
    if (censor) {
        trans.p <- censoring(data, censor.val, censor.range,
            missing.range, sigma, knn$nn.index, cb)
    }
    else {
        # パッケージ名を指定
        d2 <- switch(distance, euclidean = destiny:::d2_no_censor(knn$nn.index,
            knn$nn.dist, cb), cosine = icos2_no_censor(knn$nn.index,
            imputed.data, cb))
        # パッケージ名を指定
        trans.p <- Matrix:::sparseMatrix(d2@i, p = d2@p, x = exp(-d2@x/(2 *
            sigma^2)), dims = dim(d2), index1 = FALSE)
        rm(d2)
    }
    if (verbose) {
        close(pb)
        cat("...done. Time:\n")
        print(proc.time() - tic)
    }
    rm(knn)
    diag(trans.p) <- 0
    # パッケージ名を指定
    trans.p <- Matrix:::drop0(trans.p)
    trans.p <- Matrix:::symmpart(trans.p)
    # apply
    d <- apply(trans.p, 1, sum)

    max.dist <- max(trans.p@x, na.rm = TRUE)

    # パッケージ名を指定
    destiny:::stopifsmall(max.dist)
    if (density.norm) {
        trans.p <- as(trans.p, "dgTMatrix")
        # パッケージ名を指定
        H <- Matrix:::sparseMatrix(trans.p@i, trans.p@j, x = trans.p@x/(d[trans.p@i +
            1] * d[trans.p@j + 1]), dims = dim(trans.p), index1 = FALSE)
    }
    else {
        H <- trans.p
    }
    rm(trans.p)
    # パッケージ名を指定、apply
    D.rot <- Matrix:::Diagonal(x = apply(H, 1, sum)^-0.5)
    M <- D.rot %*% H %*% D.rot
    rm(H)
    if (!is.null(.debug.env)) {
        assign("M", M, .debug.env)
        assign("D.rot", D.rot, .debug.env)
    }
    if (verbose) {
        cat("performing eigen decomposition...")
        tic <- proc.time()
    }
    # パッケージ名を指定
    eig.M <- destiny:::eig.decomp(M, n, n.eigs, TRUE)
    if (verbose) {
        cat("...done. Time:\n")
        print(proc.time() - tic)
    }
    # as.matrixを加えた
    eig.vec <- as.matrix(t(t(eig.M$vectors) %*% as.matrix(D.rot)))
    colnames(eig.vec) <- paste0("DC", seq(0, n.eigs))
    # クラスの仕様の確認が働いてめんどくさいから、リストで返すようにした
    list(eigenvalues = eig.M$values[-1], eigenvectors = eig.vec[,
        -1], sigmas = sigmas, data.env = data.env, eigenvec0 = eig.vec[,
        1], d = d, k = k, density.norm = density.norm, distance = distance,
        censor.val = censor.val, censor.range = censor.range,
        missing.range = missing.range, M = M)
}

# カテゴリカルなカーネル関数
CatKernel <- function(label, type=c("two", "one_vs_rest", "each")){
    # データ数
    N <- length(label)

    # 各クラスの集計
    tl <- table(label)

    # 集計値ベクトル
    sum.label <- sapply(label, function(x){tl[which(labels(tl)$label == x)]})

    # 2-class
    if(length(tl) == 2){
        if((type=="two")){
            ########### Equation (9) ############
            sapply(label, function(x){
                    out <- rep(0, length=length(label))
                    # 同じクラスなら 1/m_+ - 1/m
                    out[which(x == label)] <- 1 / tl[which(labels(tl)$label == x)] - 1 / N
                    # 違うクラスなら - 1/m
                    out[which(x != label)] <- - 1 / N
                    out
                })
            #####################################
        }else{
            warning("Wrong type paramter!")
        }
    }
    # multi-class
    else if(length(tl) > 2){
        if(type=="one_vs_rest"){
            ########### Equation (11) ###########
            sapply(label, function(x){
                    out <- rep(0, length=length(label))
                    # 同じクラスなら（C_i == C_j） 1/m_i
                    out[which(x == label)] <- 1 / tl[which(labels(tl)$label == x)]

                    # 違うクラスなら（C_i != C_j） 1/ (m_i - N)
                    out[which(x != label)] <- 1 / (sum.label[which(x != label)] - N)
                    out
                })
            #####################################
        }else if(type=="each"){
            ########### Equation (12) ###########
            sapply(label, function(x){
                    out <- rep(0, length=length(label))
                    # 同じクラスなら（C_i == C_j） 1/m_i
                    out[which(x == label)] <- 1 / sqrt(tl[which(labels(tl)$label == x)])

                    # 違うクラスなら（C_i != C_j） 1/ (m_i - N)
                    out[which(x != label)] <- 0
                    out
                })
            #####################################
        }else{
            warning("Wrong type parameter!")
        }
    }else{
        warning("Please confirm your label vector!")
    }
}

# ヒルベルト-シュミット独立性基準（K行列とL行列の独立性の指標）
HSIC <- function(K, L, N){
    H <- matrix(rep(-1/N), nrow=N, ncol=N)
    diag(H) <- 1 - 1/N
    sum(diag(K %*% H %*% L %*% H)) / (N-1)^2
}

# HSICを利用した特徴量抽出
FUCHIKOMA <- function(data, mode=c("Supervised", "Unsupervised"), Comp=FALSE, label=FALSE, type=FALSE, n.eigs=10){
    # データ数
    N <- ncol(data)
    # HSIC値の格納先
    HSICs <- c()
    # 削除した遺伝子
    Rejected <- c()

    # 毎回1個の遺伝子を削除
    for(i in 1:nrow(data)){
        print(i)
        # このステップでの仮のHSICs
        tmp_HSICs <- c()
        #️ このステップで見る遺伝子
        *** <-
        # データ側のグラム行列
        dif <- custom.DiffusionMap(as.ExpressionSet(data[***,]), n.eigs=n.eigs)



        # HSICを計算
        tmp_HSICs <- c(tmp_HSICs, HSIC(dif$M, L, N))
    }

    # BAHSICの最大値を格納
    HSICs <- c(HSICs, ???)
    # 削除した遺伝子を登録
    Rejected <- c(Rejected, ???)

    # 結果を出力
    list(
        DEGs = DEGs,
        Pval = Pval,
        Qval = Qval,
        BAHSICs = BAHSICs
    )
}


####################### 試しに実行 #####################
# 標準正規分布
CellA <- data.frame(matrix(rnorm(150*20), nrow=150, ncol=20))
CellB <- data.frame(matrix(rnorm(150*20), nrow=150, ncol=20))
CellC <- data.frame(matrix(rnorm(150*20), nrow=150, ncol=20))

# DEGを指定（** 要最適な作り方調査 **）
#
testdata <- data.frame(CellA, CellB, CellC)

colnames(testdata) <- c(paste0("CellA_", 1:20), paste0("CellB_", 1:20), paste0("CellC_", 1:20))
rownames(testdata) <- paste0("Gene", 1:150)

# ラベル
label <- c(rep(1, 20), rep(2, 20), rep(3, 20))

# 普通のPCA
res.pca <- prcomp(testdata)
pairs(res.pca$rotation[, 1:20], col=label)

# オブジェクト化
testdata.obj <- as.ExpressionSet(testdata)

######## Diffusion Mapを実行（第10主成分までを見る） ####
dif1 <- DiffusionMap(testdata.obj, n.eigs=10)
pairs(eigenvectors(dif1), col=label)

# 改造Ver Diffusion Map
dif2 <- custom.DiffusionMap(testdata.obj, n.eigs=10)
pairs(dif2$eigenvectors, col=label)

################### FUCHIKOMA実行 ###################
# 教師あり（クラスラベルを与える）
result1 <- FUCHIKOMA(data=testdata, mode="Supervised", label=label, type="each", n.eigs=10)

# 教師なし（指定した主成分を使う）
result2 <- FUCHIKOMA(data=testdata, mode="Unsupervised", Comp=c(1,2), n.eigs=10)









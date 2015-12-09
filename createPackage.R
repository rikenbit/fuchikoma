##############################
###### Object definition #####
##############################

#######################################################
####### グラム行列を出力するように改造したDiffusionMap ######
######### 正規化を外したほうがうまくいく可能性あり ###########
#######################################################
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
    # kd_tree, cover_tree, CR, brute
    knn <- FNN:::get.knn(imputed.data, k, algorithm = "kd_tree") # <-
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




###########################################################
################# カテゴリカルなカーネル関数 ##################
###########################################################
CatKernel <- function(label, type=c("two", "one_vs_rest", "each", "simple")){
    # データ数
    N <- length(label)

    # 各クラスの集計
    tl <- table(label)

    # 集計値ベクトル
    sum.label <- sapply(label, function(x){tl[which(labels(tl)$label == x)]})

    # any-class
    if(type == "simple"){
        ######## 化合物-タンパク質 より #########
        sapply(label, function(x){
                out <- rep(0, length=length(label))
                # 同じクラスなら（C_i == C_j） 1 / m_i
                out[which(x == label)] <- 1

                # 違うクラスなら（C_i != C_j） 1 / (m_i - N)
                out[which(x != label)] <- -1
                out
            })
        #####################################
    }else{
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
            if(type == "one_vs_rest"){
                ########### Equation (11) ###########
                sapply(label, function(x){
                        out <- rep(0, length=length(label))
                        # 同じクラスなら（C_i == C_j） 1 / m_i
                        out[which(x == label)] <- 1 / tl[which(labels(tl)$label == x)]

                        # 違うクラスなら（C_i != C_j） 1 / (m_i - N)
                        out[which(x != label)] <- 1 / (sum.label[which(x != label)] - N)
                        out
                    })
                #####################################
            }else if(type == "each"){
                ########### Equation (12) ###########
                sapply(label, function(x){
                        out <- rep(0, length=length(label))
                        # 同じクラスなら（C_i == C_j） 1 / m_i
                        out[which(x == label)] <- 1 / sqrt(tl[which(labels(tl)$label == x)])

                        # 違うクラスなら（C_i != C_j） 1 / (m_i - N)
                        out[which(x != label)] <- 0
                        out
                    })
                #####################################
            }else{
                warning("Wrong type parameter!")
            }
        }else{
            warning("Confirm your label vector!")
        }
    }
}




###########################################################
### ヒルベルト-シュミット独立性基準（K行列とL行列の独立性の指標）####
###########################################################
HSIC <- function(K, L, p.value=FALSE){
    if(!all(c(dim(K), dim(L)) == dim(K)[1])){
        stop("Inappropriate matrices are specified!\nPlease confirm the number of rows and columns.")
    }

    N <- dim(K)[1]
    H <- matrix(rep(-1/N), nrow=N, ncol=N)
    diag(H) <- 1 - 1/N
    K <- as.matrix(K)
    L <- as.matrix(L)
    # The value of HSIC
    hsic.value <- sum(diag(K %*% H %*% L %*% H)) / (N-1)^2

    # 例外処理1
    if(hsic.value < 0){
        hsic.value <- 0
    }

    # Pvalue of HSIC
    if(p.value){
        # HSICの期待値
        u_x2 <- 1/(N*(N-1)) * sum(K[upper.tri(K)])
        u_y2 <- 1/(N*(N-1)) * sum(L[upper.tri(L)])
        E_HSIC <- 1 / N * (1 + u_x2 * u_y2 - u_x2 - u_y2)

        # HSICの分散
        H <- matrix(rep(-1/N), nrow=N, ncol=N)
        diag(H) <- 1 - 1/N
        B <- H %*% K %*% H * H %*% L %*% H
        B <- B * B
        var_HSIC <- sum(2*(N-4)*(N-5)/(N*(N-1)*(N-2)*(N-3)) * t(rep(1, length=N)) %*% (B-diag(B)))

        # 例外処理2
        if(E_HSIC <= 0){
            E_HSIC <- 0.0001
        }
        if(var_HSIC <= 0){
            var_HSIC <- 0.0001
        }
        # Shape parameter（ > 0）
        Alpha <- E_HSIC^2 / var_HSIC
        # Scale paramter（ > 0）
        Beta <- N * var_HSIC / E_HSIC
        # この値がガンマ分布に従う（ > 0）
        x <- N * hsic.value

        # p値
        p_HSIC <- pgamma(x, shape=Alpha, scale=Beta, lower.tail=FALSE)
    }else{
        p_HSIC <- NA
    }
    # return
    list(HSIC=hsic.value, Pval=p_HSIC)
}




###########################################################
############# FUCHIKOMA : HSICを利用した特徴量抽出 ###########
###########################################################
FUCHIKOMA <- function(data, mode=c("Supervised", "Unsupervised"), Comp=FALSE, label=FALSE, cat.type=FALSE, n.eigs=10, algorithm=c("brute", "song"), per.rej=10, threshold=0.01){

    # 並列化準備
    registerDoParallel(detectCores())

    ############ ラベル側のグラム行列（一回のみ） ##############
    if((mode == "Supervised") && (is.vector(label))){
        L <- CatKernel(label, type=cat.type)
    }else if(mode == "Unsupervised"){
        if(is.vector(Comp)){
            EigenVecs <- custom.DiffusionMap(as.ExpressionSet(as.data.frame(t(data))), n.eigs=n.eigs)$eigenvectors[, Comp]
            L <- EigenVecs %*% t(EigenVecs)
        }else{
            warning("Specify Comp!")
        }
    }else{
        warning("Wrong mode!")
    }
    ############ ラベル側のグラム行列（終わり）#################

        # HSIC値の格納先
        HSICs <- 0
        # 削除した遺伝子の場所
        RejPosition <- c()

        ##️############### BAHSICの計算ステップ ################
        while(length(SurvPosition) > 5){
            #️ このステップで見る遺伝子（生き残り）
            SurvPosition <- setdiff(1:nrow(data), RejPosition)
            cat(paste0("### No. of remaining gene is ", length(SurvPosition), " ###\n"))

            # 共通で使うsigmaの計算
            sigma =  try(destiny::optimal.sigma(find.sigmas(as.ExpressionSet(as.data.frame(t(data[SurvPosition,]))), verbose = FALSE)), silent = TRUE)

            # 生き残り遺伝子内でのHSIC計算
            tmp_HSICs <- foreach(j = 1:length(SurvPosition), .export=c("SurvPosition", "custom.DiffusionMap", "n.eigs", "HSIC", "data", "L", "sigma"), .combine = "c") %dopar% {
                # データ側のグラム行列
                dif <- try(custom.DiffusionMap(as.ExpressionSet(as.data.frame(t(data[SurvPosition[setdiff(1:length(SurvPosition), j)],]))), n.eigs=n.eigs, sigma=sigma))
                if ('try-error' %in% class(dif)){
                    return(0)
                }else{
                    K <- dif$M
                    # HSICを計算
                    tmp_HSIC <- HSIC(K, L)$HSIC
                    return(ifelse(is.nan(tmp_HSIC), 0, tmp_HSIC))
                }
            }
            names(tmp_HSICs) <- rownames(data)[SurvPosition]

            ############### 各ステップでの最後の処理 #############
            if(algorithm == "brute"){
                # 今回一番HSICが大きくなった遺伝子
                NoRej <- 1
            }else if(algorithm == "song"){
                # 今回一番HSICが大きくなった遺伝子
                NoRej <- round(N * per.rej / 100)
            }else{
                stop("algorithm parameter is wrong!")
            }

            # 今回HSICが最大な上位n個の遺伝子
            tmp_MaxHSICs <- rev(sort(tmp_HSICs))[1:NoRej]

            # HSICsがこれまでのHSICsの最大値よりも小さくなったら打ち切り
            if(max(HSICs) - max(tmp_MaxHSIC) < threshold){
                # BAHSICの最大値を格納
                HSICs <- c(HSICs, tmp_MaxHSICs)
                # 削除した遺伝子を登録
                RejPosition <- c(RejPosition, unlist(sapply(names(tmp_MaxHSICs), function(x){which(x == rownames(data))})))
            ########### 各ステップでの最後の処理（終わり）##########
            }else{
                break
            }
        }
        ##️############ BAHSICの計算ステップ（終わり） #############

        # HSICs整形
        HSICs <- HSICs[setdiff(2:length(HSICs), which(is.na(HSICs)))]

        # 一番HSICが大きかった遺伝子以降 => DEGs
        DEGs <- HSICs[which(max(HSICs) == HSICs):length(HSICs)]

        # non-DEGs
        nonDEGs <- HSICs[1:(which(max(HSICs) == HSICs) - 1)]

        # 最後に残った遺伝子でp値を算出（Pval用）
        pval_HSICs <- foreach(j = 1:length(DEGs), .export=c("DEGs", "HSICs", "custom.DiffusionMap", "n.eigs", "HSIC", "data", "L", "sigma"), .combine = "c") %dopar% {
            # データ側のグラム行列
            dif <- try(
                custom.DiffusionMap(
                    as.ExpressionSet(
                        as.data.frame(
                            t(
                                data[
                                names(setdiff(HSICs, DEGs)),
                                ]
                                )
                            )
                        )
                    , n.eigs=n.eigs, sigma=sigma
                    )
                )
            if ('try-error' %in% class(dif)){
                return(0)
            }else{
                K <- dif$M
                # HSICのp値を計算
                HSIC(K, L, p.value=TRUE)$Pval
            }
        }

        # 結果を出力
        list(
            DEGs.HSICs = DEGs,
            DEGs.Pvals = pval_HSICs,
            All.HSICs = HSICs,
            Rej.order = rank(nonDEGs)
        )
}

##############################
######## Objects list ########
##############################
objectslist <- c(
	"custom.DiffusionMap",
	"CatKernel",
	"HSIC",
	"FUCHIKOMA"
)

##############################
########## Packaging #########
##############################
package.skeleton(name = "FUCHIKOMA", objectslist, path = ".")

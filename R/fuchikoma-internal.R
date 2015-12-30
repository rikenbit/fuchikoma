.custom.DiffusionMap <-
function (data, sigma = NULL, k = find.dm.k(nrow(data) - 1L), 
    n.eigs = min(20L, nrow(data) - 2L), density.norm = TRUE, 
    ..., distance = c("euclidean", "cosine"), censor.val = NULL, 
    censor.range = NULL, missing.range = NULL, vars = NULL, verbose = !is.null(censor.range), 
    .debug.env = NULL) 
{
    distance <- match.arg(distance)
    data.env <- new.env(parent = .GlobalEnv)
    data.env$data <- data
    data <- destiny:::extract.doublematrix(data, vars)
    imputed.data <- data
    if (any(is.na(imputed.data))) 
        imputed.data <- as.matrix(hotdeck(data, imp_var = FALSE))
    n <- nrow(imputed.data)
    if (n <= n.eigs + 1L) 
        stop(sprintf("Eigen decomposition not possible if n <U+2264> n.eigs+1 (And %s <U+2264> %s)", 
            n, n.eigs + 1L))
    if (is.null(k) || is.na(k)) 
        k <- n - 1L
    if (k >= nrow(imputed.data)) 
        stop(sprintf("k has to be < nrow(data) (And %s <U+2265> nrow(data))", 
            k))
    censor <- destiny:::test.censoring(censor.val, censor.range, 
        data, missing.range)
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
    knn <- get.knn(imputed.data, k, algorithm = "kd_tree")
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
        d2 <- switch(distance, euclidean = destiny:::d2_no_censor(knn$nn.index, 
            knn$nn.dist, cb), cosine = icos2_no_censor(knn$nn.index, 
            imputed.data, cb))
        trans.p <- sparseMatrix(d2@i, p = d2@p, x = exp(-d2@x/(2 * 
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
    trans.p <- drop0(trans.p)
    trans.p <- symmpart(trans.p)
    d <- apply(trans.p, 1, sum)
    max.dist <- max(trans.p@x, na.rm = TRUE)
    destiny:::stopifsmall(max.dist)
    if (density.norm) {
        trans.p <- as(trans.p, "dgTMatrix")
        H <- sparseMatrix(trans.p@i, trans.p@j, x = trans.p@x/(d[trans.p@i + 
            1] * d[trans.p@j + 1]), dims = dim(trans.p), index1 = FALSE)
    }
    else {
        H <- trans.p
    }
    rm(trans.p)
    D.rot <- Diagonal(x = apply(H, 1, sum)^-0.5)
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
    eig.M <- destiny:::eig.decomp(M, n, n.eigs, TRUE)
    if (verbose) {
        cat("...done. Time:\n")
        print(proc.time() - tic)
    }
    eig.vec <- as.matrix(t(t(eig.M$vectors) %*% as.matrix(D.rot)))
    colnames(eig.vec) <- paste0("DC", seq(0, n.eigs))
    list(eigenvalues = eig.M$values[-1], eigenvectors = eig.vec[, 
        -1], sigmas = sigmas, data.env = data.env, eigenvec0 = eig.vec[, 
        1], d = d, k = k, density.norm = density.norm, distance = distance, 
        censor.val = censor.val, censor.range = censor.range, 
        missing.range = missing.range, M = M)
}
.Shrink.HSIC <-
function (K, L, H, N, HSIC) 
{
    KL <- 1/N * sum(diag(K %*% H) * diag(L %*% H))
    (1 - (KL - HSIC)/((N - 2) * HSIC + KL/N)^2) * HSIC
}
.uni.fuchikoma <-
function (data, mode = c("Supervised", "Unsupervised"), Comp = NULL, 
    label = FALSE, cat.type = c("simple", "one_vs_rest", "each", 
        "two"), kernel = vanilladot(), n.eigs = 10) 
{
    mode <- match.arg(mode, c("Supervised", "Unsupervised"))
    if (!is.null(Comp) && (Comp > n.eigs)) {
        warning("Inappropriate Comp parameter!")
    }
    cat.type <- match.arg(cat.type, c("simple", "one_vs_rest", 
        "each", "two"))
    if ((n.eigs > nrow(data)) || (0 > n.eigs)) {
        warning("Inappropriate n.eigs parameter!")
    }
    if ((mode == "Supervised") && (is.vector(label))) {
        L <- CatKernel(label, type = cat.type)
    }
    else if (mode == "Unsupervised") {
        if (is.vector(Comp)) {
            DCs <- .custom.DiffusionMap(as.ExpressionSet(as.data.frame(t(data))), 
                n.eigs = n.eigs)$eigenvectors[, Comp]
            EigenVals <- .custom.DiffusionMap(as.ExpressionSet(as.data.frame(t(data))), 
                n.eigs = n.eigs)$eigenvalues[Comp]
            DCs_e <- t(apply(DCs, 1, function(x) {
                x * EigenVals
            }))
            L <- DCs_e %*% t(DCs_e)
        }
        else {
            warning("Specify Comp!")
        }
    }
    else {
        warning("Wrong mode!")
    }
    HSICs <- apply(data, 1, function(x) {
        HSIC(kernelMatrix(kernel, t(t(x))), L)
    })
    list(All.HSICs = sapply(HSICs, function(x) {
        x$HSIC
    }), All.Pvals = sapply(HSICs, function(x) {
        x$Pval
    }))
}

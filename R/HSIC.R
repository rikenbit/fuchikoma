HSIC <-
function (K, L, shrink = FALSE, type = c("gamma", "permutation"), 
    n.perm = 100) 
{
    if (!all(c(dim(K), dim(L)) == dim(K)[1])) {
        stop("\nInappropriate matrices are specified!\nPlease confirm the number of rows and columns.")
    }
    N <- dim(K)[1]
    H <- matrix(rep(-1/N), nrow = N, ncol = N)
    diag(H) <- 1 - 1/N
    K <- as.matrix(K)
    L <- as.matrix(L)
    hsic.value <- sum(diag(K %*% H %*% L %*% H))/N^2
    if (shrink) {
        hsic.value <- .Shrink.HSIC(K, L, H, N, hsic.value)
    }
    if (hsic.value < 0) {
        hsic.value <- 0
    }
    if (type == "gamma") {
        u_x2 <- 1/(N * (N - 1)) * sum(K[upper.tri(K)])
        u_y2 <- 1/(N * (N - 1)) * sum(L[upper.tri(L)])
        E_HSIC <- 1/N * (1 + u_x2 * u_y2 - u_x2 - u_y2)
        H <- matrix(rep(-1/N), nrow = N, ncol = N)
        diag(H) <- 1 - 1/N
        B <- H %*% K %*% H * H %*% L %*% H
        B <- B * B
        var_HSIC <- sum(2 * (N - 4) * (N - 5)/(N * (N - 1) * 
            (N - 2) * (N - 3)) * t(rep(1, length = N)) %*% (B - 
            diag(B)))
        if (E_HSIC <= 0) {
            E_HSIC <- 0.1
        }
        if (var_HSIC <= 0) {
            var_HSIC <- 0.1
        }
        Alpha <- E_HSIC^2/var_HSIC
        Beta <- N * var_HSIC/E_HSIC
        x <- N * hsic.value
        p_HSIC <- pgamma(x, shape = Alpha, scale = Beta, lower.tail = FALSE)
    }
    else if (type == "permutation") {
        registerDoParallel(detectCores())
        HSICs_rand <- foreach(j = 1:n.perm, .export = c("N", 
            "H", "K", "L", "shrink", ".Shrink.HSIC"), .combine = "c") %dopar% 
            {
                rand_order <- sample(1:N, N)
                L_rand <- L[rand_order, rand_order]
                hsic.value.rand <- sum(diag(K %*% H %*% L_rand %*% 
                  H))/N^2
                if (shrink) {
                  hsic.value.rand <- .Shrink.HSIC(K, L, H, N, 
                    hsic.value.rand)
                }
                if (hsic.value.rand < 0) {
                  0
                }
                else {
                  hsic.value.rand
                }
            }
        p_HSIC <- length(which(HSICs_rand > hsic.value))/n.perm
    }
    else {
        warning("Wrong type parameter!")
    }
    list(HSIC = hsic.value, Pval = p_HSIC)
}

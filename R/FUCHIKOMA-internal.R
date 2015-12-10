.Shrink.HSIC <-
function (K, L, H, N, HSIC) 
{
    KL <- 1/N * sum(diag(K %*% H) * diag(L %*% H))
    (1 - (KL - HSIC)/((N - 2) * HSIC + KL/N)^2) * HSIC
}

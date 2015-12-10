.Shrink.HSIC <-
function (K, L, H, N, HSIC) 
{
    KL <- 1/N * sum(diag(K %*% H) * diag(L %*% H))
    (1 - (KL - hsic.value)/((N - 2) * hsic.value + KL/N)^2) * 
        hsic.value
}

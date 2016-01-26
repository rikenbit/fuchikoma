simulateDropoutCounts <-
function (Ngene = 10000, PDEG = 0.2, DEG.assign = c(0.5, 0.5), 
    DEG.foldchange = c(4, 4), replicates = c(3, 3), Lambda = 0.1) 
{
    if (!is.numeric(Ngene)) {
        warning("Please specify the Ngene as numeric!")
    }
    if ((PDEG > 1) || (0 > PDEG)) {
        warning("Please specify the PDEG as numeric (0 < PDEG < 1)!")
    }
    if (sum(DEG.assign) != 1) {
        warning("Please specify the DEG.assign as vector(sum(DEG.assign) is 1)!")
    }
    if ((length(DEG.assign) != length(DEG.foldchange)) || (length(DEG.foldchange) != 
        length(replicates))) {
        warning("Please specify the DEG.assign, DEG.foldchange, and DEG.foldchange as same length vectors!")
    }
    if (!is.numeric(Lambda)) {
        warning("Please specify the Lambda as numeric!")
    }
    Nreplicate <- sum(replicates)
    No.DEG <- Ngene * PDEG * DEG.assign
    data(Marioni)
    mu <- apply(Marioni, 1, mean)
    v <- apply(Marioni, 1, var)
    Disp <- (v - mu)/mu^2
    rn.index <- sample(which(Disp > 0), Ngene, replace = TRUE)
    design.matrix <- matrix(1, nrow = Ngene, ncol = length(replicates))
    for (i in 1:length(replicates)) {
        if (i == 1) {
            row.index <- 1:No.DEG[1]
        }
        else {
            row.index <- sum(No.DEG[1:(i - 1)]) + 1:sum(No.DEG[i])
        }
        design.matrix[row.index, i] <- DEG.foldchange[i]
    }
    original.matrix <- t(sapply(rn.index, function(x) {
        rnbinom(n = Nreplicate, mu = mu[x], size = 1/Disp[x])
    }))
    for (i in 1:length(replicates)) {
        deg.index <- which(design.matrix[, i] != 1)
        col.index <- sum(replicates[i - 1]) + 1:sum(replicates[i])
        original.matrix[deg.index, col.index] <- t(sapply(rn.index[deg.index], 
            function(x) {
                rnbinom(n = replicates[i], mu = DEG.foldchange[i] * 
                  mu[x], size = 1/Disp[x])
            }))
    }
    mean.vector <- apply(original.matrix, 1, mean)
    var.vector <- apply(original.matrix, 1, var)
    droprate <- exp(-Lambda * mean.vector^2)
    droprate.matrix <- sapply(1:Nreplicate, function(y) {
        sapply(droprate, function(x) {
            rbinom(1, 1, prob = (1 - x))
        })
    })
    testdata.matrix <- original.matrix * droprate.matrix
    colnames(testdata.matrix) <- paste0("Cell_", 1:ncol(testdata.matrix))
    rownames(testdata.matrix) <- paste0("Gene_", 1:nrow(testdata.matrix))
    return(list(DEG = rownames(testdata.matrix)[which(rowMeans(design.matrix) != 
        1)], FC = design.matrix, Simcount = testdata.matrix))
}

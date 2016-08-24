simulateDropoutCounts <-
function (Ngene = 10000, makeDEG = TRUE, PDEG = 0.02, DEG.assign = c(0.5, 
    0.5), DEG.thr = c("E5", "E5"), replicates = c(3, 3), Lambda = 0.1) 
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
    if ((length(DEG.assign) != length(DEG.thr))) {
        warning("Please specify the DEG.assign and DEG.thr are the same length!")
    }
    if (!is.numeric(Lambda)) {
        warning("Please specify the Lambda as numeric!")
    }
    data(Marioni)
    mu <- apply(Marioni, 1, mean)
    v <- apply(Marioni, 1, var)
    Disp <- (v - mu)/mu^2
    rn.index <- sample(which(Disp > 0), Ngene, replace = TRUE)
    original.matrix <- t(sapply(rn.index, function(x) {
        rnbinom(n = sum(replicates), mu = mu[x], size = 1/Disp[x])
    }))
    if (makeDEG) {
        No.DEG <- Ngene * PDEG * DEG.assign
        design.matrix <- matrix(1, nrow = Ngene, ncol = length(replicates))
        for (i in 1:length(replicates)) {
            if (i == 1) {
                row.index <- 1:No.DEG[1]
            }
            else {
                row.index <- sum(No.DEG[1:(i - 1)]) + 1:sum(No.DEG[i])
            }
            design.matrix[row.index, i] <- .GenerateFC(mu[rn.index][row.index], 
                DEG.thr[i])
        }
        for (i in 1:length(replicates)) {
            deg.index <- which(design.matrix[, i] != 1)
            col.index <- sum(replicates[1:i - 1]) + 1:sum(replicates[i])
            original.matrix[deg.index, col.index] <- t(sapply(deg.index, 
                function(x) {
                  rnbinom(n = replicates[i], mu = design.matrix[x, 
                    i] * mu[rn.index[x]], size = 1/Disp[rn.index[x]])
                }))
        }
    }
    mean.vector <- apply(original.matrix, 1, mean)
    var.vector <- apply(original.matrix, 1, var)
    droprate <- exp(-Lambda * mean.vector^2)
    droprate.matrix <- sapply(1:sum(replicates), function(y) {
        sapply(droprate, function(x) {
            rbinom(1, 1, prob = (1 - x))
        })
    })
    testdata.matrix <- original.matrix * droprate.matrix
    colnames(testdata.matrix) <- as.vector(unlist(sapply(1:length(replicates), 
        function(x) {
            paste0("Group", x, "_rep", 1:replicates[x])
        })))
    rownames(testdata.matrix) <- paste0("Gene_", 1:nrow(testdata.matrix))
    Group <- unlist(sapply(1:length(replicates), function(x) {
        group <- rep(x, length = replicates[x])
        names(group) <- rep(paste0("Group", x), length = replicates[x])
        group
    }, simplify = FALSE))
    if (makeDEG) {
        return(list(Group = Group, DEG = rownames(testdata.matrix)[which(rowMeans(design.matrix) != 
            1)], FC = design.matrix, Simcount = testdata.matrix))
    }
    else {
        return(list(Group = Group, Simcount = testdata.matrix))
    }
}

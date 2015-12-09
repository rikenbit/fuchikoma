FUCHIKOMA <-
function (data, mode = c("Supervised", "Unsupervised"), Comp = FALSE, 
    label = FALSE, cat.type = FALSE, n.eigs = 10, algorithm = c("song", 
        "brute"), per.rej = 10, threshold = 0.01) 
{
    registerDoParallel(detectCores())
    if ((mode == "Supervised") && (is.vector(label))) {
        L <- CatKernel(label, type = cat.type)
    }
    else if (mode == "Unsupervised") {
        if (is.vector(Comp)) {
            EigenVecs <- custom.DiffusionMap(as.ExpressionSet(as.data.frame(t(data))), 
                n.eigs = n.eigs)$eigenvectors[, Comp]
            L <- EigenVecs %*% t(EigenVecs)
        }
        else {
            warning("Specify Comp!")
        }
    }
    else {
        warning("Wrong mode!")
    }
    HSICs <- 0
    All.pval <- 0
    RejPosition <- c()
    SurvPosition <- 1:nrow(data)
    while (length(SurvPosition) > 5) {
        cat(paste0("### No. of remaining gene is ", length(SurvPosition), 
            " ###\n"))
        sigma = try(destiny::optimal.sigma(find.sigmas(as.ExpressionSet(as.data.frame(t(data[SurvPosition, 
            ]))), verbose = FALSE)), silent = TRUE)
        if ("try-error" %in% class(sigma)) {
            cat(paste0("Error in destiny::optimal.sigma !!\n"))
            break
        }
        tmp_HSICs_Pvals <- foreach(j = 1:length(SurvPosition), 
            .export = c("SurvPosition", "custom.DiffusionMap", 
                "n.eigs", "HSIC", "data", "L", "sigma")) %dopar% 
            {
                dif <- try(custom.DiffusionMap(as.ExpressionSet(as.data.frame(t(data[SurvPosition[setdiff(1:length(SurvPosition), 
                  j)], ]))), n.eigs = n.eigs, sigma = sigma))
                if ("try-error" %in% class(dif)) {
                  return(NA)
                }
                else {
                  K <- dif$M
                  HSIC(K, L, p.value = TRUE)
                }
            }
        tmp_HSICs <- unlist(lapply(tmp_HSICs_Pvals, function(x) {
            x$HSIC
        }))
        tmp_Pvals <- unlist(lapply(tmp_HSICs_Pvals, function(x) {
            x$Pval
        }))
        names(tmp_HSICs) <- rownames(data)[SurvPosition]
        names(tmp_Pvals) <- rownames(data)[SurvPosition]
        if (algorithm == "brute") {
            NoRej <- 1
        }
        else if (algorithm == "song") {
            NoRej <- round(length(SurvPosition) * per.rej/100)
        }
        else {
            stop("algorithm parameter is wrong!")
        }
        tmp_MaxHSICs <- rev(sort(tmp_HSICs))[1:NoRej]
        if (max(HSICs) - max(tmp_MaxHSICs) < threshold) {
            HSICs <- c(HSICs, tmp_MaxHSICs)
            All.pval <- c(All.pval, tmp_Pvals[names(tmp_MaxHSICs)])
            RejPosition <- c(RejPosition, unlist(sapply(names(tmp_MaxHSICs), 
                function(x) {
                  which(x == rownames(data))
                })))
            SurvPosition <- setdiff(1:nrow(data), RejPosition)
        }
        else {
            break
        }
    }
    HSICs <- HSICs[setdiff(2:length(HSICs), which(is.na(HSICs)))]
    All.pval <- All.pval[setdiff(2:length(All.pval), which(is.na(All.pval)))]
    DEGs <- HSICs[which(max(HSICs) == HSICs):length(HSICs)]
    nonDEGs <- HSICs[1:(which(max(HSICs) == HSICs) - 1)]
    list(DEGs.HSICs = DEGs, DEGs.Pvals = All.pval[names(DEGs)], 
        All.HSICs = HSICs, All.Pvals = All.pval, Rej.order = sort(rank(nonDEGs)))
}

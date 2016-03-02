fuchikoma <-
function (data, cores = NULL, mode = c("Supervised", "Unsupervised", 
    "Mix"), weight = c(0.5, 0.5), Comp = NULL, label = FALSE, 
    cat.type = c("simple", "one_vs_rest", "each", "two"), n.eigs = 10, 
    algorithm = c("song", "brute"), per.rej = 10, threshold = 0.01, 
    verbose = FALSE, dropout = 10, sigma = 60) 
{
    if (!is.null(cores) && (cores < 1)) {
        warning("Inappropriate cores parameter!")
    }
    mode <- match.arg(mode, c("Supervised", "Unsupervised", "Mix"))
    if (sum(weight) != 1) {
        warning("Sum of weight must be 1!")
    }
    if (!is.null(Comp) && (Comp > n.eigs)) {
        warning("Inappropriate Comp parameter!")
    }
    cat.type <- match.arg(cat.type, c("simple", "one_vs_rest", 
        "each", "two"))
    if ((n.eigs > nrow(data)) || (0 > n.eigs)) {
        warning("Inappropriate n.eigs parameter!")
    }
    algorithm <- match.arg(algorithm, c("song", "brute"))
    if ((per.rej > 100) || (0 > per.rej)) {
        warning("Inappropriate per.rej parameter!")
    }
    if (!is.numeric(threshold)) {
        warning("Inappropriate threshold parameter!")
    }
    if ((dropout > 100) || (0 > dropout)) {
        warning("Inappropriate dropout parameter!")
    }
    if (is.null(cores)) {
        registerDoParallel(detectCores())
    }
    else {
        registerDoParallel(cores)
    }
    on.exit(stopImplicitCluster())
    L <- .Lmatrix(data, mode = mode, weight = weight, Comp = Comp, 
        label = label, cat.type = cat.type, n.eigs = n.eigs, 
        sigma = sigma)
    HSICs <- 0
    All.pval <- 0
    RejPosition <- c()
    SurvPosition <- 1:nrow(data)
    counter <- 0
    while (length(SurvPosition) > dropout) {
        if (verbose) {
            counter <- counter + 1
            try(dir.create("verbose"))
            try(dir.create(paste0("verbose/step_", counter)))
            cat(paste0("### No. of remaining gene is ", length(SurvPosition), 
                " ###\n"))
        }
        tmp_HSICs_Pvals <- foreach(j = 1:length(SurvPosition), 
            .export = c("SurvPosition", ".custom.DiffusionMap", 
                "n.eigs", "HSIC", "data", "L", "sigma", "verbose", 
                "counter")) %dopar% {
            dif <- try(.custom.DiffusionMap(as.ExpressionSet(as.data.frame(t(data[SurvPosition[setdiff(1:length(SurvPosition), 
                j)], ]))), n.eigs = n.eigs, sigma = sigma))
            if ("try-error" %in% class(dif)) {
                return(NA)
            }
            else {
                K <- dif$M
                HSIC(K, L)
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
        if (verbose) {
            print("###### HSIC ######")
            print(head(tmp_HSICs))
            print("###### P-value ######")
            print(head(tmp_Pvals))
        }
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
    remaining <- names(rev(sort(tmp_HSICs))[(NoRej + 1):length(tmp_HSICs)])
    tmp_HSICs <- tmp_HSICs[remaining]
    tmp_Pvals <- tmp_Pvals[remaining]
    HSICs <- c(HSICs[setdiff(2:length(HSICs), which(is.na(HSICs)))], 
        tmp_HSICs)
    All.pval <- c(All.pval[setdiff(2:length(All.pval), which(is.na(All.pval)))], 
        tmp_Pvals)
    if (max(HSICs) != mean(HSICs)) {
        DEGs <- HSICs[which(max(HSICs) == HSICs):length(HSICs)]
    }
    else {
        DEGs <- NA
    }
    list(DEGs.HSICs = DEGs, DEGs.Pvals = All.pval[names(DEGs)], 
        All.HSICs = HSICs, All.Pvals = All.pval)
}

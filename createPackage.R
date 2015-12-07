################################# Definition of Objects ####################################
load("BreastCancer.rda")
BreastCancer <- data
load("pvals.rda")
pvals <- output
load("StudyA.rda")
load("result.rda")
Result.Meta <- result

############################################################################################

original.busca <- function (x, S){
	which(S[,1] == x[1] & S[,2] == x[2])
}
original.n.menor <- function (x, S1, S2){
	length(which(S1 <= x[1] &  S2 <= x[2]))
}

text.busca_unix <- "
# include <Rcpp.h>
// [[Rcpp::export]]

Rcpp::NumericVector busca (Rcpp::NumericVector x, Rcpp::NumericMatrix S){
using std::vector;
vector<int> out;
int rS = S.nrow();

for(int i=0; i<rS; i++){
	if((S(i,0) == x[0]) && (S(i,1) == x[1])){
		int j = i + 1;
		out.push_back(j);
		break;
	}
}
return Rcpp::wrap(out);
}
"

text.n.menor_unix <- "
# include <Rcpp.h>

// [[Rcpp::export]]

Rcpp::NumericVector nmenor (Rcpp::NumericVector x, Rcpp::NumericVector S1, Rcpp::NumericVector S2){

using std::vector;
vector<int> out;
int rS = S1.size();

for(int i=0; i<rS; i++){
	if((S1[i] <= x[0]) && (S2[i] <= x[1])){
		out.push_back(i);
	}
}
return Rcpp::wrap(out.size());
}
"

# Windowsにも対応できるように、改行文字を書き換え
text.n.menor_win <- gsub("\n", "\r\n", text.n.menor_unix)
text.busca_win <- gsub("\n", "\r\n", text.busca_unix)

Accelerate.NOISeq <- function(OS=NULL){
	if(OS == "Unix"){
		data(text.n.menor_unix)
		data(text.busca_unix)
		env <- getNamespace("NOISeq")
		sourceCpp(code = text.n.menor_unix)
		sourceCpp(code = text.busca_unix)
		assignInNamespace("busca", busca, ns="NOISeq", envir=env)
		assignInNamespace("n.menor", nmenor, ns="NOISeq", envir=env)		
	}
	if(OS == "Windows"){
		data(text.n.menor_win)
		data(text.busca_win)
		env <- getNamespace("NOISeq")
		sourceCpp(code = text.n.menor_win)
		sourceCpp(code = text.busca_win)
		assignInNamespace("busca", busca, ns="NOISeq", envir=env)
		assignInNamespace("n.menor", nmenor, ns="NOISeq", envir=env)
	}
	else{
		warning("Please specify OS as Unix or Windows!")
	}
}

Reset.Accelerate.NOISeq <- function(){
	env <- getNamespace("NOISeq")
	assignInNamespace("busca", original.busca, ns="NOISeq", envir=env)
	assignInNamespace("n.menor", original.n.menor, ns="NOISeq", envir=env)
}

############################################################################################

meta.readData <- function(data = NULL, factors = NULL, length = NULL, biotype = NULL, chromosome = NULL, gc = NULL, studies = NULL){
	if(is.null(studies)){
		stop("Please specify \"studies\" at first!\n")
	}
	else if(!is.vector(factors)){
		stop("Please \"factors\" parameter as vector.\n")
	}
	else if(!is.vector(studies)){
		stop("Please \"studies\" parameter as vector.\n")
	}
	else if(length(factors) != length(studies)){
		stop("Length of factors and that of studies are different!\n")
	}else{

		out <- list()
		l <- nlevels(as.factor(studies))
		length(out) <- l

		length2 <- NULL
		biotype2 <- NULL
		chromosome2 <- NULL
		gc2 <- NULL

		if(!is.null(length)){length2 <- length}
		if(!is.null(biotype)){biotype2 <- biotype}
		if(!is.null(chromosome)){chromosome2 <- chromosome}
		if(!is.null(gc)){gc2 <- gc}

		e <<- new.env()
		e$factors <- factors
		e$studies <- studies
		e$length <- length2
		e$biotype <- biotype2
		e$chromosome <- chromosome2
		e$gc <- gc2

		loc <- list()
		length(loc) <- l
		e$loc <- loc

		for(x in 1:l){
			loc[[x]] <- which(e$studies == levels(as.factor(e$studies))[x])
			e$loc[[x]] <- loc[[x]]
			out[[x]] <- readData(data = data[, e$loc[[x]]], factors = as.data.frame(e$factors[e$loc[[x]]]),
				length = e$length, biotype = e$biotype, chromosome = e$chromosome, gc = e$gc)
		}
	}
	class(out) <- "metaExpressionSet"
	return(out)
}

############################################################################################

original.probdeg <- function (Mg, Dg, Mn, Dn, prec = 2) {
    tot <- length(Mn)
    gens <- names(Mg)
    Mruido <- abs(round(Mn, prec))
    Druido <- round(Dn, prec)
    Mgen <- abs(round(Mg, prec))
    Dgen <- round(Dg, prec)
    MDgen <- na.omit(cbind(Mgen, Dgen))
    MDunic <- unique(MDgen)
    Nres <- apply(MDunic, 1, n.menor, S1 = Mruido, S2 = Druido)
    lugares <- apply(MDgen, 1, busca, S = MDunic)
    Nconj <- Nres[lugares]
    names(Nconj) <- names(lugares)
    laprob <- Nconj/tot
    laprob <- laprob[gens]
    names(laprob) <- gens
    Nconj <- Nconj[gens]
    names(Nconj) <- gens
    laprob <- list(prob = laprob, numDE = Nconj, numNOISE = tot)
    laprob
}

############################################################################################

original.MD <- function (dat = dat, selec = c(1:nrow(dat))) {
    pares <- as.matrix(combn(ncol(dat), 2))
    if (NCOL(pares) > 30) {
        sub30 <- sample(1:NCOL(pares), size = 30, replace = FALSE)
        pares <- pares[, sub30]
    }
    mm <- NULL
    dd <- NULL
    for (i in 1:ncol(pares)) {
        a <- dat[selec, pares[1, i]]
        b <- dat[selec, pares[2, i]]
        mm <- cbind(mm, log(a/b, 2))
        dd <- cbind(dd, abs(a - b))
    }
    list(M = mm, D = dd)
}
############################################################################################

custom.probdeg <- function (Mg, Dg, Mn, Dn, prec = 2) {
	tot <- length(Mn)
	gens <- names(Mg)
	Mruido <- round(Mn, prec)
	Druido <- round(Dn, prec)
	Mgen <- round(Mg, prec)
	Dgen <- round(Dg, prec)
	MDgen <- na.omit(cbind(Mgen, Dgen))
	MDunic <- unique(MDgen)  
	Nres <- apply(MDunic, 1, NOISeq:::n.menor, S1 = Mruido, S2 = Druido)
	lugares <- apply(MDgen, 1, NOISeq:::busca, S = MDunic)
	Nconj <- Nres[lugares]
	names(Nconj) <- names(lugares)
	laprob <- Nconj / tot
	laprob <- laprob[gens]
	names(laprob) <- gens
	Nconj <- Nconj[gens]
	names(Nconj) <- gens
	laprob <- list("prob" = laprob, "numDE" = Nconj, "numNOISE" = tot)
	laprob
}

############################################################################################

custom.MD <- function (dat = dat, selec = c(1:nrow(dat))) {
  pares <- as.matrix(combn(ncol(dat), 2))

  if (NCOL(pares) > 30) {  
    sub30 <- sample(1:NCOL(pares), size = 30, replace = FALSE)
    pares <- pares[,sub30]
  }
  
  mm <- NULL
  dd <- NULL
  for (i in 1:ncol(pares)) {
    a <- dat[selec,pares[1,i]]
    b <- dat[selec,pares[2,i]]
    mm <- cbind(mm, log(a/b, 2))
    dd <- cbind(dd, (a-b))
  }
  list("M" = mm, "D" = dd)
}

############################################################################################

oneside.noiseq <- function(input, k = 0.5, norm = c("rpkm", "uqua", "tmm", "n"), 
    replicates = c("technical", "biological", "no"), factor = NULL, 
    conditions = NULL, pnr = 0.2, nss = 5, v = 0.02, lc = 1, x = NULL){

	env <- getNamespace("NOISeq")
	assignInNamespace("probdeg", custom.probdeg, ns="NOISeq", envir=env)
	assignInNamespace("MD", custom.MD, ns="NOISeq", envir=env)

	k2 = 0.5
	norm2 = c("rpkm", "uqua", "tmm", "n")
   	replicates2 = c("technical", "biological", "no")
   	factor2 = NULL
   	conditions2 = NULL
   	pnr2 = 0.2
   	nss2 = 5
   	v2 = 0.02
   	lc2 = 1

   	if(k != 0.5){k2 <- k}
   	if(length(norm) != 4){norm2 <- norm}
   	if(length(replicates) != 3){replicates2 <- replicates}
   	if(!is.null(factor)){factor2 <- factor}
   	if(!is.null(conditions)){conditions2 <- conditions}
   	if(pnr != 0.2){pnr2 <- pnr}
   	if(nss != 5){nss2 <- nss}
   	if(v != 0.02){v2 <- v}
   	if(lc != 1){lc2 <- lc}

  	e2 <<- new.env()
  	e2$k <- k2
  	e2$norm <- norm2
  	e2$replicates <- replicates2
  	e2$factor <- factor2
  	e2$conditions <- conditions2
  	e2$pnr <- pnr2
  	e2$nss <- nss2
  	e2$v <- v2
  	e2$lc <- lc2

	# Return
	fff <- eval(parse(text=eval(parse(text="e2$factor"))))
	A <- length(which(conditions[1] == fff))
	B <- length(which(conditions[2] == fff))

	if((A==1)&&(B==1)){
		warning("Your dataset contains some non-replicated sample. \"replicates\" parameter is automatically selected as \"no\".\n")
		out <- NOISeq::noiseq(input, k = e2$k, norm = e2$norm, 
    	replicates = "no", factor = e2$factor, 
    	conditions = e2$conditions, pnr = e2$pnr, nss = e2$nss, v = e2$v, lc = e2$lc)
	}else{
	out <- NOISeq::noiseq(input, k = e2$k, norm = e2$norm, 
    replicates = e2$replicates, factor = e2$factor, 
    conditions = e2$conditions, pnr = e2$pnr, nss = e2$nss, v = e2$v, lc = e2$lc)
	}
	return(out)

	# Unload
	assignInNamespace("probdeg", original.probdeg, ns="NOISeq", envir=env)
	assignInNamespace("MD", original.MD, ns="NOISeq", envir=env)
}

############################################################################################

meta.oneside.noiseq <- function(input, k = 0.5, norm = c("rpkm","uqua","tmm","n"),
			replicates = c("technical","biological","no"),
			factor = NULL, conditions = NULL, pnr = 0.2, nss = 5, v = 0.02, lc = 1, studies = NULL, cl = NULL){
	
	if (inherits(input, "metaExpressionSet") == FALSE) {
        stop("Error. You must give an object generated by the meta.readData function or other.oneside.pvalues\n")
	}
	if(is.null(studies)){
		stop("Please specify \"studies\" at first!\n")
	}
	else if(!is.vector(factor)){
		stop("Please \"factor\" parameter as vector.\n")
	}
	else if(!is.vector(studies)){
		stop("Please \"studies\" parameter as vector.\n")
	}
	else if(length(factor) != length(studies)){
		stop("Length of factors and that of studies are different!\n")
	}else{
		out <- list()
		l <- length(input)
		length(out) <- l

		k2 = 0.5
		norm2 = c("rpkm", "uqua", "tmm", "n")
	   	replicates2 = c("technical", "biological", "no")
	   	conditions2 = NULL
	   	pnr2 = 0.2
	   	nss2 = 5
	   	v2 = 0.02
	   	lc2 = 1

	   	if(k != 0.5){k2 <- k}
	   	if(length(norm) != 4){norm2 <- norm}
	   	if(length(replicates) != 3){replicates2 <- replicates}
	   	if(!is.null(conditions)){conditions2 <- conditions}
	   	if(pnr != 0.2){pnr2 <- pnr}
	   	if(nss != 5){nss2 <- nss}
	   	if(v != 0.02){v2 <- v}
	   	if(lc != 1){lc2 <- lc}

	  	e <<- new.env()
	  	e$input <- input
		e$factors <- factor
		e$studies <- studies
		e$k <- k2
		e$norm <- norm2
		e$replicates <- replicates2
		e$conditions <- conditions2
		e$pnr <- pnr2
		e$nss <- nss2
		e$v <- v2
		e$lc <- lc2

		loc <- list()
		length(loc) <- l
		e$loc <- loc	
		for(x in 1:l){
			e$loc[[x]] <- which(levels(as.factor(e$studies))[x] == e$studies)	
		}

		# Return
		if(!is.null(cl)){
			clusterExport(cl, "e")
			#clusterEvalQ(cl, library("NOISeq"))
			out <- snow::parSapply(cl, 1:l, function(x){
					output <- list()
					length(output) <- 3
					names(output) <- c("upper", "lower", "weight")
					result <- metaSeq:::oneside.noiseq(e$input[[x]], k = e$k, norm = e$norm, replicates = e$replicates, factor = "e$factors[e$loc[[x]]]", conditions = e$conditions, pnr = e$pnr, nss = e$nss, v = e$v, lc = e$lc, x = x)
					U <- result@results[[1]]$prob
					L <- 1 - result@results[[1]]$prob
					W <- nrow(input[[x]]@phenoData@data)
					names(U) <- rownames(input[[x]]@assayData$exprs)
					names(L) <- rownames(input[[x]]@assayData$exprs)
					output$upper <- U
					output$lower <- L
					output$weight <- W
					return(output)
					}
				   )

			colnames(out) <- paste("Study", 1:l)
			return(out)	
		}else{
			out <- sapply(1:l, function(x){
					output <- list()
					length(output) <- 3
					names(output) <- c("upper", "lower", "weight")
					result <- metaSeq:::oneside.noiseq(e$input[[x]], k = e$k, norm = e$norm, replicates = e$replicates, factor = "e$factors[e$loc[[x]]]", conditions = e$conditions, pnr = e$pnr, nss = e$nss, v = e$v, lc = e$lc, x = x)
					U <- result@results[[1]]$prob
					L <- 1 - result@results[[1]]$prob
					W <- nrow(input[[x]]@phenoData@data)
					names(U) <- rownames(input[[x]]@assayData$exprs)
					names(L) <- rownames(input[[x]]@assayData$exprs)
					output$upper <- U
					output$lower <- L
					output$weight <- W
					return(output)
					}
				   )

			colnames(out) <- paste("Study", 1:l)
			return(out)	
		}
	}
}

############################################################################################

Fisher.test <- function(pvals, na.mode = "notignore"){
	l <- ncol(pvals)
	U <- pvals[1,1][[1]]
	L <- pvals[2,1][[1]]
	weight <- unlist(pvals["weight",])

	for(i in 2:l){
		U <- cbind(U, pvals[1,i][[1]])
		L <- cbind(L, pvals[2,i][[1]])
	}

	up <- c()
	low <- c()

	if(na.mode == "notignore"){
		up <- apply(U, 1, each.Fisher.test)
		low <- apply(L, 1, each.Fisher.test)
	}
	if(na.mode == "ignore"){
		up <- apply(U, 1, each.Fisher.ignore.test)
		low <- apply(L, 1, each.Fisher.ignore.test)
	}else{
		warnings("You have to specify na.mode as \"ignore\" or \"notignore\"")
	}
	list("Upper" = up, "Lower" = low, "Weight" = weight)
}

############################################################################################

each.Fisher.test <- function(p) {
  Xsq <- -2*sum(log(p))
  pval <- 1-pchisq(Xsq, df = 2*length(p))
  return(pval)
}

############################################################################################

each.Fisher.ignore.test <- function(p) {
  loc <- setdiff(1:length(p), which(is.na(p)))
  if(length(loc) >= 1){
	  p <- p[loc]
  	Xsq <- -2*sum(log(p))
  	pval <- 1-pchisq(Xsq, df = 2*length(p))
  	return(pval)
  }else{
  	return(NA)
  }
}

############################################################################################

Stouffer.test <- function(pvals, na.mode = "notignore"){

	if(is.null(pvals["weight",][[1]])){
		stop("Weight is needed for Stouffer's method!\n")
	}else{
		l <- ncol(pvals)
		U <- pvals[1,1][[1]]
		L <- pvals[2,1][[1]]
		weight <- unlist(pvals["weight",])

		for(i in 2:l){
			U <- cbind(U, pvals[1,i][[1]])
			L <- cbind(L, pvals[2,i][[1]])
		}

		up <- c()
		low <- c()

		if(na.mode == "notignore"){
			up <- apply(U, 1, each.Stouffer.test, w=weight)
			low <- apply(L, 1, each.Stouffer.test, w=weight)
		}
		else if(na.mode == "ignore"){
			up <- apply(U, 1, each.Stouffer.ignore.test, w=weight)
			low <- apply(L, 1, each.Stouffer.ignore.test, w=weight)
		}else{
			warnings("You have to specify na.mode as \"ignore\" or \"notignore\"")
		}
		list("Upper" = up, "Lower" = low, "Weight" = weight)
	}
}

############################################################################################

each.Stouffer.test <- function(p, w) {
  if (missing(w)) {
    w <- rep(1, length(p))/length(p)
  } else {
    if (length(w) != length(p))
      stop("Length of p and w must equal!")
  }
  Zi <- qnorm(1-p) 
  Z  <- sum(w*Zi)/sqrt(sum(w^2))
  pval <- 1-pnorm(Z)
  return(pval)
}

############################################################################################

each.Stouffer.ignore.test <- function(p, w) {
  if (missing(w)) {
    w <- rep(1, length(p))/length(p)
  } else {
    if (length(w) != length(p))
      stop("Length of p and w must equal!")
  }

  loc <- setdiff(1:length(p), which(is.na(p)))
  if(length(loc) >= 1){
	p <- p[loc]
	w <- w[loc]
	Zi <- qnorm(1-p) 
  	Z  <- sum(w*Zi)/sqrt(sum(w^2))
  	pval <- 1-pnorm(Z)
  	return(pval)
  }else{
  	return(NA)
  }
}

############################################################################################

other.oneside.pvalues <- function (Upper, Lower, weight = NULL) 
{
    if ((min(!is.na(Upper)) < 0) || (max(!is.na(Upper)) > 1)) {
        stop("Is this dataset (upper) pvalues? Some elements exceed 0 - 1 range. Please confirm first.\n")
    }
    if ((min(!is.na(Lower)) < 0) || (max(!is.na(Lower)) > 1)) {
        stop("Is this dataset (lower) pvalues? Some elements exceed 0 - 1 range. Please confirm first.\n")
    }
    if (nrow(Upper) != nrow(Lower)) {
        stop("Number of rows in upper p-values and lower p-values are different!\n")
    }
    if (ncol(Upper) != ncol(Lower)) {
        stop("Number of columns in upper p-values and lower p-values are different!\n")
    }
    if ((!is.null(weight)) && (ncol(Upper) != length(weight))) {
        A <- ncol(Upper)
        B <- length(weight)
        Call <- paste0("Number of column in p-value matrix is ", 
            A, " but length of weight vector is ", B, ". Please confirm first.\n")
        stop(Call)
    }
    l <- ncol(Upper)
    out <- sapply(1:l, function(x) {
        output <- list()
        length(output) <- 3
        names(output) <- c("upper", "lower", "weight")
		U <- Upper[, x]
		L <- Lower[, x]
		names(U) <- rownames(Upper)
		names(L) <- rownames(Lower)
		output$upper <- U
		output$lower <- L
		output$weight <- weight[x]
		return(output)
    })
    colnames(out) <- paste("Exp", 1:l)
    return(out)
}

############################################################################################


# Objects list
objectslist <- c(
	"BreastCancer",
	"pvals",
	"StudyA",
	"Result.Meta",
	"meta.readData",
	"original.probdeg",
	"original.MD",
	"custom.probdeg",
	"custom.MD",
	"oneside.noiseq",
	"meta.oneside.noiseq",
	"Fisher.test",
	"each.Fisher.test",
	"each.Fisher.ignore.test",
	"Stouffer.test",
	"each.Stouffer.test",
	"each.Stouffer.ignore.test",
	"other.oneside.pvalues",
	# 2013/11/16 追加分	
	"original.n.menor",
	"original.busca",
	# "text.n.menor",
	# "text.busca",
	"Accelerate.NOISeq",
	"Reset.Accelerate.NOISeq",
	# 2013/12/2 追加分	
	"text.n.menor_unix",
	"text.busca_unix",
	"text.n.menor_win",
	"text.busca_win"
)

# Packaging
package.skeleton(name = "metaSeq", objectslist, path = ".")

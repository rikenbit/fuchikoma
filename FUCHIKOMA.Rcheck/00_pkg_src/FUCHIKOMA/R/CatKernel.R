CatKernel <-
function (label, type = c("two", "one_vs_rest", "each", "simple")) 
{
    N <- length(label)
    tl <- table(label)
    sum.label <- sapply(label, function(x) {
        tl[which(labels(tl)$label == x)]
    })
    if (type == "simple") {
        sapply(label, function(x) {
            out <- rep(0, length = length(label))
            out[which(x == label)] <- 1
            out[which(x != label)] <- -1
            out
        })
    }
    else {
        if (length(tl) == 2) {
            if ((type == "two")) {
                sapply(label, function(x) {
                  out <- rep(0, length = length(label))
                  out[which(x == label)] <- 1/tl[which(labels(tl)$label == 
                    x)] - 1/N
                  out[which(x != label)] <- -1/N
                  out
                })
            }
            else {
                warning("Wrong type paramter!")
            }
        }
        else if (length(tl) > 2) {
            if (type == "one_vs_rest") {
                sapply(label, function(x) {
                  out <- rep(0, length = length(label))
                  out[which(x == label)] <- 1/tl[which(labels(tl)$label == 
                    x)]
                  out[which(x != label)] <- 1/(sum.label[which(x != 
                    label)] - N)
                  out
                })
            }
            else if (type == "each") {
                sapply(label, function(x) {
                  out <- rep(0, length = length(label))
                  out[which(x == label)] <- 1/sqrt(tl[which(labels(tl)$label == 
                    x)])
                  out[which(x != label)] <- 0
                  out
                })
            }
            else {
                warning("Wrong type parameter!")
            }
        }
        else {
            warning("Confirm your label vector!")
        }
    }
}

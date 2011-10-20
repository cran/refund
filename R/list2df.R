list2df <-
function(l){
    nrows <- sapply(l, function(x) nrow(as.matrix(x)))
    stopifnot(length(unique(nrows)) == 1)
    ret <- data.frame(rep(NA, nrows[1]))
    for(i in 1:length(l)) ret[[i]] <- l[[i]]
    names(ret) <- names(l)
    return(ret)
}


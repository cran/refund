summary.pffr <-
function (object, ...) {
    call <- match.call()
    call[[1]] <- mgcv::summary.gam
    class(object) <- class(object)[-1]
    call$object <- object
    ret <- eval(call)
    
    ret$formula <- object$pffr$formula
    
    shrtlbls <- getShrtlbls(object)
    
    rownames(ret$s.table) <- sapply(rownames(ret$s.table), 
            function(x){
                shrtlbls[pmatch(x, unlist(object$pffr$labelmap))]     
            })
    class(ret) <- c("summary.pffr", class(ret))
    ret$n  <- paste(ret$n, "(", object$pffr$nobs," x ", object$pffr$nyindex, ")", sep="") 
    return(ret)
}


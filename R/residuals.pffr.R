residuals.pffr <-
function (object, reformat=TRUE, ...) 
{
    if (!inherits(object, "pffr")) 
        stop("`object' is not of class \"pffr\"")
    ret <- mgcv:::residuals.gam(object, ...)
    if(reformat) ret <- matrix(ret, nrow=object$pffr$nobs, ncol=object$pffr$nyindex, byrow=TRUE)
    return(ret)
}


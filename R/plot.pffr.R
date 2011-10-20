plot.pffr <- function (x, ...)
{
    call <- match.call()
    call[[1]] <- mgcv::plot.gam
    #drop "pffr" class and replace <x> with changed value s.t. method dispatch works without glitches
    class(x) <- class(x)[-1]
    call$x <- x
    invisible(eval(call))
}


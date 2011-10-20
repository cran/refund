model.matrix.pffr <-
function (object, ...) 
{
    if (!inherits(object, "pffr")) 
        stop("`object' is not of class \"pffr\"")
    predict(object, type = "lpmatrix", reformat=FALSE, ...)
}


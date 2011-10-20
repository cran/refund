expand.call <-
function(definition=NULL, call=sys.call(sys.parent(1)), expand.dots = TRUE)
{
	call <- .Internal(match.call(definition, call, expand.dots))
	ans <- as.list(call)
	
	frmls <- formals(safeDeparse(ans[[1]]))
	frmls <- frmls[!sapply(frmls, is.symbol)]
	
	add <- which(!(names(frmls) %in% names(ans)))
	return(as.call(c(ans, frmls[add])))
}


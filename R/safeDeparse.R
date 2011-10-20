safeDeparse <-
function(expr){
	ret <- paste(deparse(expr), collapse="")
	gsub("[[:space:]][[:space:]]+", " ", ret)
}


getShrtlbls <-
function(object){
    
    labelmap <- object$pffr$labelmap
    
    ret <- sapply(names(unlist(labelmap)),
            function(x){
                x <- gsub("((?:[A-Za-z]+))(\\s)(=)(\\s)", "", x, perl=TRUE)
                x <- gsub("\\s", "", x)
                if(length(gregexpr(",",x)[[1]])>=3) {
                    x <- gsub("([^,]*,[^,]*,[^,]*)(.*)(\\)$)", "\\1\\3", x, perl=TRUE)    
                }
                if(any(regexpr(",\".*",x)[[1]]>0)) {
                    x <- gsub("([^\"]*)(,[c\\(]*\".*)(\\)$)", "\\1\\3", x, perl=TRUE) 
                }
                return(x)
            })
    if(any(sapply(labelmap, length) > 1)){
        which <- which(sapply(labelmap, length) > 1)
        inds <- c(0, cumsum(sapply(labelmap, length)))
        for(w in which){
            ret[(inds[w]+1):inds[w+1]] <- {
                lbls <- labelmap[[w]]
                bylevels <- sapply(object$smooth[lbls], function(x) x$by.level)
                by <- object$smooth[[lbls[1]]]$by
                paste(by, bylevels, "(", object$pffr$yindname, ")", sep="")
            }
        }
    }        
    if(any(!grepl("\\(", ret))){
        which <- which(!grepl("\\(", ret))
        ret[which] <- paste(ret[which],"(", object$pffr$yindname, ")", sep="")
    }
    
    return(ret)
}


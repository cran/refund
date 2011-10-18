predict.pffr <-
function(object,
        newdata,
        reformat=TRUE,
        type = "link",
        se.fit = FALSE,
        terms = NULL, 
        ...){
    
    call <- match.call()
    nyindex <- object$pffr$nyindex
    
    if(!is.null(terms)) {
        warning("<terms> argument not implemented yet, ignored.")
        terms <- NULL
    }
    
    if(!missing(newdata)){
        if(!(all(names(newdata) %in% names(object$model)))){
            stopifnot(length(unique(sapply(newdata, function(x) ifelse(is.matrix(x), nrow(x), length(x))))) ==1)
            
            nobs <- nrow(as.matrix(newdata[[1]]))
            
            gamdata <- list()
            gamdata[[paste(object$pffr$yindname, ".vec", sep="")]] <- rep(object$pffr$yind, times=nobs)
            
            varmap <- sapply(names(object$pffr$labelmap), function(x) all.vars(formula(paste("~", x))))
            
            for(cov in names(newdata)){
                trms <- which(sapply(varmap, function(x) any(grep(paste("^",cov,"$",sep=""), x)))) 
                if(length(trms)!=0){
                    for(trm in trms){
                        is.ff <- trm %in% object$pffr$where$where.ff
                        is.sff <- trm %in% object$pffr$where$where.sff
                        if(is.ff){
                            ff <- object$pffr$ff[[grep(paste(cov,"[,\\)]",sep=""), names(object$pffr$ff))]]
                            if(grepl(paste(cov,"\\.[st]mat",sep=""), deparse(ff$call$x))){
                                L <- ff$L 
                                if(any(apply(L, 2, function(x) length(unique(x)))!=1)){
                                    stop("Error for ", names(varmap)[trm],
                                            "-- Prediction for ff-terms with varying rows in integration operator L not implememented yet.")
                                } 
                                predL <- matrix(L[1,], byrow=TRUE, nrow=nrow(newdata[[cov]]), ncol=ncol(L))
                                
                                gamdata[[paste(cov, ".smat", sep="")]] <- 
                                        matrix(ff$xind, byrow=TRUE, ncol=length(ff$xind), nrow=nobs*nyindex)
                                gamdata[[paste(cov, ".tmat", sep="")]] <- 
                                        matrix(rep(ff$yind, times=nobs), ncol=length(ff$xind), nrow=nobs*nyindex)
                                gamdata[[paste("L.", cov, sep="")]] <-  
                                        (predL*newdata[[cov]])[rep(1:nobs, e=nyindex),]
                            }
                        }
                        if(is.sff){
                            sff <- object$pffr$ff[[grep(paste("^",cov,"$",sep=""), names(object$pffr$ff))]]
                            if(grepl(paste(cov,"\\.[st]mat",sep=""), deparse(sff$call$x))){
                                L <- sff$L 
                                if(any(apply(L, 2, function(x) length(unique(x)))!=1)){
                                    stop("Error for ", names(varmap)[trm],
                                            "-- Prediction for sff-terms with varying rows in integration operator L not implememented yet.")
                                } 
                                predL <- matrix(L[1,], byrow=TRUE, nrow=nrow(newdata[[cov]]), ncol=ncol(L))
                                
                                gamdata[[paste(cov, ".mat", sep="")]] <- newdata[[cov]][rep(1:nobs, e=nyindex),]
                                gamdata[[paste(cov, ".smat", sep="")]] <- 
                                        matrix(sff$xind, byrow=TRUE, ncol=length(sff$xind), nrow=nobs*nyindex)
                                gamdata[[paste(cov, ".tmat", sep="")]] <- 
                                        matrix(rep(sff$yind, times=nobs), ncol=length(sff$xind), nrow=nobs*nyindex)
                                gamdata[[paste("L.", cov, sep="")]] <-  predL[rep(1:nobs, e=nyindex),]
                            }
                        }
                        if(!(is.ff || is.sff)) {
                            gamdata[[cov]] <- drop(newdata[[cov]])[rep(1:nobs, e=nyindex)]
                        }
                        
                    } 
                }
            }
            gamdata <- list2df(gamdata)
            call[["newdata"]] <- gamdata
        }
    } else {
        call$newdata <- eval(call$newdata)
        nobs <- object$pffr$nobs  
    }
    
    
    call[[1]] <- mgcv::predict.gam
    ret <- eval(call)
    
    if(type=="lpmatrix" && reformat){
        reformat <- FALSE
        warning("Setting reformat to FALSE for type=\"lpmatrix\".")
    }
    
    if(reformat){
        if(se.fit){
            if(type %in% c("terms", "iterms")){
                ret <- lapply(ret, function(x)
                            do.call(list,
                                    sapply(1:ncol(x), function(i){
                                                d <- list(I(matrix(x[,i], nrow=nobs, ncol=object$pffr$nyindex, byrow=TRUE)))
                                                names(d)  <- colnames(x)[i]
                                                return(d)
                                            })))
            } else {
                ret <- lapply(ret, function(x) matrix(x, nrow=nobs, ncol=object$pffr$nyindex, byrow=TRUE))
            } 
        } else {
            if(type %in% c("terms", "iterms")){
                ret <- do.call(list, sapply(1:ncol(ret), function(i){
                                    d <- list(I(matrix(ret[,i], nrow=nobs, ncol=object$pffr$nyindex, byrow=TRUE)))
                                    names(d)  <- colnames(ret)[i]
                                    return(d)
                                }))
            } else ret <- matrix(ret, nrow=nobs, ncol=object$pffr$nyindex, byrow=TRUE)
        }
        if(!missing(terms) && (type %in% c("terms", "iterms"))){
            ind <- sapply(terms, function(x)
                        sapply(x, function(x)
                                    which(x == unlist(object$pffr$labelmap))))
            attr(ret, "names") <- names(object$pffr$labelmap)[ind]
        }
    }
    return(ret)
}


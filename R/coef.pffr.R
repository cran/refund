coef.pffr <- function(object, raw=FALSE, se=TRUE, freq=FALSE, n1=100, n2=40, n3=20, ...){
    if(raw){
        return(object$coefficients)
    } else {
        getCoefs <- function(i){
            ## this constructs a grid over the range of the covariates
            ## and returns estimated values on this grid, with
            ## by-variables set to 1
            ## cf. mgcv:::plots.R (plot.mgcv.smooth etc..) for original code
            safeRange <- function(x){
                if(is.factor(x)) return(c(NA, NA))
                return(range(x, na.rm=TRUE))
            }


            makeDataGrid <- function(trm){
                #generate grid of values in range of original data
                if(trm$dim==1){
                    x <- get.var(trm$term, object$model)
                    xg <- if(is.factor(x)) {
                                unique(x)
                            } else seq(min(x), max(x), length=n1)
                    d <- data.frame(xg)
                    colnames(d) <- trm$term
                    attr(d, "xm") <- xg
                }
                if(trm$dim > 1){
                    varnms <- unlist(sapply(trm$margin, function(mar) mar$term))
                    ng <- ifelse(trm$dim==2, n2, n3)

                    x <- get.var(varnms[1], object$model)
                    xg <- if(is.factor(x)) {
                                unique(x)
                            } else seq(min(x), max(x),length=ng)
                    y <- get.var(varnms[2], object$model)
                    yg <- if(is.factor(y)) {
                                unique(y)
                            } else seq(min(y), max(y),length=ng)
                    if(length(varnms)==2){
                        d <- expand.grid(xg, yg)
                        attr(d, "xm") <- xg
                        attr(d, "ym") <- yg
                    } else {
                        z <- get.var(varnms[3], object$model)
                        zg <- if(is.factor(z)) {
                                    unique(z)
                                } else seq(min(z), max(z), length=ng)
                        d <- expand.grid(xg, yg, zg)
                        attr(d, "xm") <- xg
                        attr(d, "ym") <- yg
                        attr(d, "zm") <- zg
                    }
                    colnames(d) <- varnms

                }
                if(trm$by!="NA"){
                    d$by <- 1
                    colnames(d) <- c(head(colnames(d),-1), trm$by)
                }
                return(d)
            }


            getP <- function(trm, d){
                #return an object similar to what plot.mgcv.smooth etc. return
                X <- PredictMat(trm, d)
                P <- if(trm$dim==1){
                            list(x=attr(d, "xm"), xlab=trm$term, xlim=safeRange(attr(d, "xm")))
                        } else {
                            varnms <- unlist(sapply(trm$margin, function(mar) mar$term))
                            if(length(varnms) == 2){
                                list(x=attr(d, "xm"), y=attr(d, "ym"), xlab=varnms[1], ylab=varnms[2],
                                        ylim=safeRange(attr(d, "ym")), xlim=safeRange(attr(d, "xm")))
                            } else {
                                if(trm$dim==3){
                                    list(x=attr(d, "xm"), y=attr(d, "ym"), z=attr(d, "zm"),
                                            xlab=varnms[1], ylab=varnms[2], zlab=varnms[3],
                                            ylim=safeRange(attr(d, "ym")), xlim=safeRange(attr(d, "xm")), zlim=safeRange(attr(d, "zm")))
                                }
                            }
                        }
                trmind <- trm$first.para:trm$last.para
                P$value <- X%*%object$coefficients[trmind]
                P$coef <- cbind(d, "value"=P$value)
                if(se){
                    P$se <- sqrt(rowSums((X%*%covmat[trmind, trmind])*X))
                    P$coef <- cbind(P$coef, se=P$se)
                }
                P$dim <- trm$dim
                return(P)
            }

            trm <- object$smooth[[i]]
            if(trm$dim > 3){
                warning("can't deal with smooths with more than 3 dimensions, returning NULL for ", shrtlbls[i])
                return(NULL)
            }

            d <- makeDataGrid(trm)
            P <- getP(trm, d)

            #browser()
            # get proper labeling
            P$main <- shrtlbls[names(object$smooth)[i] == unlist(object$pffr$labelmap)]
            which <- match(names(object$smooth)[i], object$pffr$labelmap)
            if(which %in% object$pffr$where$where.ff){
                which.ff <- which(object$pffr$where$where.ff == which)
                P$xlab <- object$pffr$yindname
                ylab <- deparse(as.call(formula(paste("~",names(object$pffr$ff)[which.ff]))[[2]])$xind)
                if(ylab=="NULL") ylab <- "xindex"
                P$ylab <- ylab
            }
            if(which %in% object$pffr$where$where.sff){

                which.sff <- which(object$pffr$where$where.sff == which)
                P$xlab <- object$pffr$yindname
                ylab <- deparse(as.call(formula(paste("~",names(object$pffr$ff)[which.sff]))[[2]])$xind)
                if(ylab=="NULL") ylab <- "xindex"
                P$ylab <- ylab
                P$zlab <- gsub(".mat$", "", object$pffr$ff[[which.sff]]$xname)
            }

            return(P)
        }

        covmat <- if(freq){
                    object$Ve
                } else {
                    object$Vp
                }
        ret <- list()
        smind <- unlist(sapply(object$smooth, function(x){
                            seq(x$first.para, x$last.para)
                        }))
        ret$pterms <- cbind(value=object$coefficients[-smind])
        if(se) ret$pterms <- cbind(ret$pterms, se=sqrt(diag(covmat)[-smind]))

        shrtlbls <- getShrtlbls(object)

        ret$smterms <- lapply(1:length(object$smooth), getCoefs)
        names(ret$smterms) <- sapply(seq_along(ret$smterms), function(i){
                    ret$smterms[[i]]$main
                })
        return(ret)
    }
}


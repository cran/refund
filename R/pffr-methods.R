# methods for pffr-objects
#
#
# Author: fabians
# 16.08.2011, 13:01:24
###############################################################################

#' Prediction for penalized function-on-function regression
#'
#'  Takes a fitted \code{pffr}-object produced by \code{\link{pffr}()} and produces
#'  predictions given a new set of values for the model covariates or the original
#'  values used for the model fit. Predictions can be accompanied by standard errors,
#'  based on the posterior distribution of the model coefficients. This is a wrapper
#'  function for \code{\link[mgcv]{predict.gam}()}.
#'
#'  Index variables (i.e., evaluation points) for the functional covariates are reused
#'  from the fitted model object and cannot be supplied with \code{newdata}.
#'  Prediction is always for the entire index range of the responses as defined
#'  in the original fit. If the original fit was performed on sparse or irregular,
#'  non-gridded response data supplied via \code{pffr}'s \code{ydata}-argument
#'  and no \code{newdata} was supplied, this function will
#'  simply return fitted values for the original evaluation points of the response (in list form).
#'  If the original fit was performed on sparse or irregular data and \code{newdata} \emph{was}
#'  supplied, the function will return predictions on the grid of evaluation points given in
#'  \code{object$pffr$yind}.
#'
#' @param object a fitted \code{pffr}-object
#' @param newdata  A named list (or a \code{data.frame}) containing the values of the
#' model covariates at which predictions are required.
#' If no \code{newdata} is provided then predictions corresponding to the original data
#' are returned. If \code{newdata} is provided then it must contain all the variables needed
#' for prediction, in the format supplied to \code{pffr}, i.e., functional predictors must be
#'  supplied as matrices with each row corresponding to one observed function.
#'  See Details for more on index variables and prediction for models fit on
#'  irregular or sparse data.
#' @param reformat logical, defaults to TRUE. Should predictions be returned in matrix form (default) or
#' in the long vector shape returned by \code{predict.gam()}?
#' @param type see \code{\link[mgcv]{predict.gam}()} for details.
#'  Note that \code{type == "lpmatrix"} will force \code{reformat} to FALSE.
#' @param se.fit see \code{\link[mgcv]{predict.gam}()}
#' @param ...  additional arguments passed on to \code{\link[mgcv]{predict.gam}()}
#' @seealso \code{\link[mgcv]{predict.gam}()}
#' @return If \code{type == "lpmatrix"}, the design matrix for the supplied covariate values in long format.
#'  If \code{se == TRUE}, a list with entries \code{fit} and \code{se.fit} containing fits and standard errors, respectively.
#'  If \code{type == "terms"} or \code{"iterms"} each of these lists is a list of matrices of the same dimension as the response for \code{newdata}
#'  containing the linear predictor and its se for each term.
#' @export
#' @method predict pffr
#' @author Fabian Scheipl
#' @importFrom mgcv predict.gam predict.bam
predict.pffr <- function(object,
                         newdata,
                         reformat=TRUE,
                         type = "link",
                         se.fit = FALSE,
                         ...){
  #browser()

  call <- match.call()
  nyindex <- object$pffr$nyindex

  ## warn if any entries in ... are not arguments for predict.gam
  dots <- list(...)
  if(length(dots)){
    validDots <- c(names(formals(predict.gam)), "cluster")
    # should be
    # unique(c(names(formals(predict.gam)),
    #          names(formals(predict.bam))))
    # but predict.bam is not exported.
    notUsed <- names(dots)[!(names(dots) %in% validDots)]
    if(length(notUsed))
      warning("Arguments <", paste(notUsed, collapse=", "), "> supplied but not used." )
  }


  if(!missing(newdata)){
    nobs <- nrow(as.matrix(newdata[[1]]))

    # check if the supplied data already has the shape expected by predict.gam
    # and dispatch immediately if so (need this so summary works as expected!)
    if(!(all(names(newdata) %in% names(object$model))) |
       !(paste0(object$pffr$yindname,".vec") %in% names(newdata))){
      # check lengths
      stopifnot(length(unique(sapply(newdata, function(x)
        ifelse(is.matrix(x), nrow(x), length(x))))) ==1)
      #        #FIXME: better leave this check to predict.gam....
      #        covnames <- mapply(gsub,
      #                pattern=c(".[st]mat$"),
      #                replacement="", x=unique(unlist(sapply(object$smooth, function(x) x$term))))
      #        covnames <- unique(covnames[covnames != paste(object$pffr$yindname, ".vec", sep="")])
      #        stopifnot(all(covnames %in% names(newdata)))


      #get newdata into the shape expected by predict gam:
      gamdata <- list()
      #y-index
      gamdata[[paste(object$pffr$yindname, ".vec", sep="")]] <- rep(object$pffr$yind, times=nobs)

      # which covariates occur in which terms?
      varmap <- sapply(names(object$pffr$labelmap), function(x) all.vars(formula(paste("~", x))))

      # don't include response
      covnames <- unique(names(newdata)[names(newdata)!=deparse(object$formula[[2]])])
      for(cov in covnames){
        #find the term(s) <cov> is associated with
        trms <- which(sapply(varmap, function(x) any(grep(paste("^",cov,"$",sep=""), x))))
        if(!is.null(dots$terms)) trms <- trms[names(trms) %in% dots$terms]
        if(length(trms)!=0){
          for(trm in trms){
            is.ff <- trm %in% object$pffr$where$ff
            is.sff <- trm %in% object$pffr$where$sff
            is.ffpc <- trm %in% object$pffr$where$ffpc
            is.pcre <- trm %in% object$pffr$where$pcre
            #if ff(X) or sff(X), generate (X.mat), X.tmat, X.smat, L.X ...
            if(is.ff){
              ff <- object$pffr$ff[[grep(paste(cov,"[,\\)]",sep=""), names(object$pffr$ff))]]
              #... but don't generate new data unless <cov> is the functional covariate.
              if(grepl(paste(cov,"\\.[st]mat",sep=""), deparse(ff$call$x))){
                # make L-matrix for new obs:
                L <- ff$L
                if(any(apply(L, 2, function(x) length(unique(x)))!=1)){
                  stop("Error for ", names(varmap)[trm],
                       "-- Prediction for ff-terms with varying rows in integration operator L not implememented yet.")
                }
                if(!is.null(ff$limits)){
                  #TODO implement prediction with limits
                  stop("Error for ", names(varmap)[trm],
                       "-- Prediction for ff-terms with <limits> not implememented yet.")
                }

                predL <- matrix(L[1,], byrow=TRUE, nrow=nrow(newdata[[cov]]), ncol=ncol(L))


                gamdata[[paste(cov, ".smat", sep="")]] <-
                  matrix(ff$xind, byrow=TRUE, ncol=length(ff$xind), nrow=nobs*nyindex)
                gamdata[[paste(cov, ".tmat", sep="")]] <-
                  matrix(rep(object$pffr$yind, times=nobs), ncol=length(ff$xind), nrow=nobs*nyindex)
                gamdata[[paste("L.", cov, sep="")]] <-
                  (predL*newdata[[cov]])[rep(1:nobs, each=nyindex),]
              }
            }
            if(is.sff){
              sff <- object$pffr$ff[[grep(paste(cov,"[,\\)]",sep=""), names(object$pffr$ff))]]
              #... but don't generate new data unless <cov> is the functional covariate.
              if(grepl(paste(cov,"\\.[st]mat",sep=""), deparse(sff$call$x))){
                # make L-matrix for new obs:
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
                  matrix(rep(object$pffr$yind, times=nobs), ncol=length(sff$xind), nrow=nobs*nyindex)
                gamdata[[paste("L.", cov, sep="")]] <-  predL[rep(1:nobs, e=nyindex),]
              }
            }
            if(is.pcre){
              pcre <- object$pffr$pcre[[grep(cov, names(object$pffr$pcre))]]
              gamdata[[paste(cov, ".vec", sep="")]] <- rep(newdata[[cov]], each=nyindex)
              for(nm in colnames(pcre$efunctions)){
                tmp <- approx(x=pcre$yind,
                              y=pcre$efunctions[, nm],
                              xout=object$pffr$yind,
                              method = "linear")$y
                gamdata[[nm]] <- tmp[rep(1:nyindex, times=nobs)]
              }
            }
            if(is.ffpc){
              ffpc <- object$pffr$ffpc[[grep(paste(cov,"[,\\)]",sep=""),
                                             names(object$pffr$ffpc))]]
              # Xc' = Phi xi' + error --> get loadings for new data:
              Xct <- t(newdata[[cov]]) - as.vector(ffpc$meanX)
              xiMat <- t(qr.coef(qr(ffpc$PCMat), Xct))
              colnames(xiMat) <- paste(make.names(cov),".PC", 1:ncol(xiMat), sep="")
              xiMat <- xiMat[rep(1:nobs, each = nyindex), , drop = FALSE]
              for(nm in colnames(xiMat)){
                gamdata[[nm]] <- xiMat[,nm]
              }
            }
            if(!(is.ff | is.sff | is.ffpc | is.pcre)) {
              gamdata[[cov]] <- if (!is.matrix(drop(newdata[[cov]]))) {
                #just repeat each entry nyindex-times to correspond to vec(<Response>)
                drop(newdata[[cov]])[rep(1:nobs, each=nyindex)]
              } else {
                # stack the matrix of the functional covariate (row-wise!)
                as.vector(t(newdata[[cov]]))
              }

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
  isIrregular <- missing(newdata) & object$pffr$sparseOrNongrid


  #call predict.gam
  call[[1]] <- if(inherits(object, "bam")){
    mgcv::predict.bam
  }  else mgcv::predict.gam
  call$object <- as.name("object")
  ret <- eval(call)

  if(type=="lpmatrix" && reformat){
    reformat <- FALSE
    warning("Setting reformat to FALSE for type=\"lpmatrix\".")
  }


  #reformat into matrices with same shape as <Response>

  if(reformat){
    if(!isIrregular){
      if(missing(newdata) && !is.null(object$pffr$missingind)){
        #pad with NAs at the appropriate locations so that fits are nobs x nyindex:
        insertNA <- function(x){
          if(length(x) != nobs*object$pffr$nyindex){
            tmp <- rep(NA, nobs*object$pffr$nyindex)
            tmp[-object$pffr$missingind] <- x
            return(tmp)
          } else {
            return(x)
          }
        }
      } else insertNA <- function(x) return(x)

      if(se.fit){
        if(type %in% c("terms", "iterms")){
          ret <- lapply(ret, function(x)
            do.call(list,
                    sapply(1:ncol(x), function(i){
                      #browser()
                      d <- list(I(matrix(insertNA(x[,i]), nrow=nobs,
                                         ncol=object$pffr$nyindex,
                                         byrow=TRUE)))
                      names(d)  <- colnames(x)[i]
                      return(d)
                    })))

        } else {
          ret <- lapply(ret, function(x) matrix(insertNA(x), nrow=nobs,
                                                ncol=object$pffr$nyindex, byrow=TRUE))
        }
      } else {
        if(type %in% c("terms", "iterms")){
          ret <- do.call(list, sapply(1:ncol(ret), function(i){
            #browser()
            d <- list(I(matrix(insertNA(ret[,i]), nrow=nobs,
                               ncol=object$pffr$nyindex, byrow=TRUE)))
            names(d)  <- colnames(ret)[i]
            return(d)
          }))
        } else ret <- matrix(insertNA(ret), nrow=nobs, ncol=object$pffr$nyindex, byrow=TRUE)
      }
    } else {
      evalpoints <- object$pffr$ydata[,c(".obs", ".index")]
      if(se.fit){
        if(type %in% c("terms", "iterms")){
          ret <- lapply(ret, function(x)
            do.call(list,
                    sapply(1:ncol(x), function(i){
                      #browser()
                      d <- list(cbind(evalpoints, .value=x[,i]))
                      names(d)  <- colnames(x)[i]
                      return(d)
                    })))

        } else {
          ret <- lapply(ret, function(x) cbind(evalpoints, .value=x))
        }
      } else {
        if(type %in% c("terms", "iterms")){
          ret <- do.call(list, sapply(1:ncol(ret), function(i){
            #browser()
            d <- list(cbind(evalpoints, .value=ret[,i]))
            names(d)  <- colnames(ret)[i]
            return(d)
          }))
        } else ret <- cbind(evalpoints, .value=ret)
      }
    }
  }
  return(ret)
}

#' Obtain model matrix for a pffr fit
#'
#' @param object a fitted \code{pffr}-object
#' @param ... other arguments, passed to \code{\link[mgcv]{predict.gam}}.
#'
#' @return A model matrix
#' @method model.matrix pffr
#' @author Fabian Scheipl
model.matrix.pffr <- function (object, ...)
{
  if (!inherits(object, "pffr"))
    stop("`object' is not of class \"pffr\"")
  predict(object, type = "lpmatrix", reformat=FALSE, ...)
}

#' Obtain residuals and fitted values for a pffr models
#'
#' See \code{\link{predict.pffr}} for alternative options to extract estimated
#' values from a \code{pffr} object.
#' "Fitted values" here refers to the estimated additive predictor values,
#' these will not be on the scale of the response for models with link functions.
#'
#' @param object a fitted \code{pffr}-object
#' @param reformat logical, defaults to TRUE. Should residuals be returned in
#'   \code{n x yindex} matrix form (regular grid data) or, respectively, in the
#'   shape of the originally supplied \code{ydata} argument (sparse/irregular
#'   data), or, if \code{FALSE}, simply as a long vector as returned by
#'   \code{resid.gam()}?
#' @param ... other arguments, passed to \code{\link[mgcv]{residuals.gam}}.
#'
#' @return A matrix or \code{ydata}-like \code{data.frame} or a vector of
#'   residuals / fitted values (see \code{reformat}-argument)
#' @export
#' @importFrom mgcv residuals.gam
#' @method residuals pffr
#' @aliases fitted.pffr
#' @author Fabian Scheipl
residuals.pffr <- function (object, reformat=TRUE, ...)
{
  if (!inherits(object, "pffr"))
    stop("`object' is not of class \"pffr\"")
  ret <- mgcv::residuals.gam(object, ...)
  if(reformat){
    if(!object$pffr$sparseOrNongrid){
      if(!(length(ret)==object$pffr$nobs*object$pffr$nyindex)){
        tmp <- rep(NA, object$pffr$nobs*object$pffr$nyindex)
        tmp[-object$pffr$missingind] <- ret
        ret <- tmp
      }
      ret <- matrix(ret, nrow=object$pffr$nobs, ncol=object$pffr$nyindex, byrow=TRUE)
    } else {
      tmp <- object$pffr$ydata
      tmp[,".value"] <- ret
      ret <- tmp
    }
  }
  return(ret)
}

#' @method fitted pffr
#' @export
#' @rdname residuals.pffr
fitted.pffr <- function (object, reformat=TRUE, ...)
{
  if (!inherits(object, "pffr"))
    stop("`object' is not of class \"pffr\"")
  ret <- object$fitted.values
  if(reformat){
    if(!object$pffr$sparseOrNongrid){
      if(!(length(ret)==object$pffr$nobs*object$pffr$nyindex)){
        tmp <- rep(NA, object$pffr$nobs*object$pffr$nyindex)
        tmp[-object$pffr$missingind] <- ret
        ret <- tmp
      }
      ret <- matrix(ret, nrow=object$pffr$nobs, ncol=object$pffr$nyindex, byrow=TRUE)
    } else {
      tmp <- object$pffr$ydata
      tmp[,".value"] <- ret
      ret <- tmp
    }
  }
  return(ret)
}

#' Plot a pffr fit
#'
#' Plot a fitted pffr-object. Simply dispatches to \code{\link[mgcv]{plot.gam}}.
#'
#' @param x a fitted \code{pffr}-object
#' @param ... arguments handed over to \code{\link[mgcv]{plot.gam}}
#'
#' @return This function only generates plots.
#' @method plot pffr
#' @importFrom mgcv plot.gam
#' @author Fabian Scheipl
plot.pffr <- function (x, ...)
{
  call <- match.call()
  call[[1]] <- mgcv::plot.gam
  #drop "pffr" class and replace <object> with changed value s.t. method dispatch works without glitches
  class(x) <- class(x)[-1]
  invisible(eval(call))
}


#' Get estimated coefficients from a pffr fit
#'
#' Returns estimated coefficient functions/surfaces \eqn{\beta(t), \beta(s,t)}
#' and estimated smooth effects \eqn{f(z), f(x,z)} or \eqn{f(x, z, t)} and their point-wise estimated standard errors.
#' Not implemented for smooths in more than 3 dimensions.
#'
#' The \code{seWithMean}-option corresponds to the \code{"iterms"}-option in \code{\link[mgcv]{predict.gam}}.
#' The \code{sandwich}-options works as follows: Assuming that the residual vectors \eqn{\epsilon_i(t), i=1,\dots,n} are i.i.d.
#' realizations of a mean zero Gaussian process with covariance \eqn{K(t,t')}, we can construct an estimator for
#' \eqn{K(t,t')} from the \eqn{n} replicates of the observed residual vectors. The covariance matrix of the stacked observations
#' vec\eqn{(Y_i(t))} is then given by a block-diagonal matrix with \eqn{n} copies of the estimated \eqn{K(t,t')} on the diagonal.
#' This block-diagonal matrix is used to construct the "meat" of a sandwich covariance estimator, similar to Chen et al. (2012),
#' see reference below.
#'
#'
#' @param object a fitted \code{pffr}-object
#' @param raw logical, defaults to FALSE. If TRUE, the function simply returns \code{object$coefficients}
#' @param se logical, defaults to TRUE. Return estimated standard error of the estimates?
#' @param freq logical, defaults to FALSE. If FALSE, use posterior variance \code{object$Vp} for variability estimates,
#'  else use \code{object$Ve}. See \code{\link[mgcv]{gamObject}}
#' @param sandwich logical, defaults to FALSE. Use a Sandwich-estimator for approximate variances? See Details.
#' 	THIS IS AN EXPERIMENTAL FEATURE, USE A YOUR OWN RISK.
#' @param seWithMean logical, defaults to TRUE. Include uncertainty about the intercept/overall mean in  standard errors returned for smooth components?
#' @param n1 see below
#' @param n2 see below
#' @param n3 \code{n1, n2, n3} give the number of gridpoints for 1-/2-/3-dimensional smooth terms
#' used in the marginal equidistant grids over the range of the covariates at which the estimated effects are evaluated.
#' @param Ktt (optional) an estimate of the covariance operator of the residual process \eqn{\epsilon_i(t) \sim N(0, K(t,t'))},
#' evaluated on \code{yind} of \code{object}. If not supplied, this is estimated from the crossproduct matrices of the
#' observed residual vectors. Only relevant for sandwich CIs.
#' @param ... other arguments, not used.
#'
#' @return If \code{raw==FALSE}, a list containing \itemize{
#'  \item \code{pterms} a matrix containing the parametric / non-functional coefficients (and, optionally, their se's)
#'  \item \code{smterms} a named list with one entry for each smooth term in the model. Each entry contains
#'     \itemize{
#'          \item \code{coef} a matrix giving the grid values over the covariates, the estimated effect (and, optionally, the se's).
#'                          The first covariate varies the fastest.
#'          \item \code{x, y, z} the unique gridpoints used to evaluate the smooth/coefficient function/coefficient surface
#'          \item \code{xlim, ylim, zlim} the extent of the x/y/z-axes
#'          \item \code{xlab, ylab, zlab} the names of the covariates for the x/y/z-axes
#'          \item \code{dim} the dimensionality of the effect
#'          \item \code{main} the label of the smooth term (a short label, same as the one used in \code{summary.pffr})
#' }}
#' @references Chen, H., Wang, Y., Paik, M.C., and Choi, A. (2013).
#' A marginal approach to reduced-rank penalized spline smoothing with application to multilevel functional data.
#' \emph{Journal of the American Statistical Association}, 101, 1216--1229.
#' @method coef pffr
#' @export
#' @importFrom mgcv PredictMat get.var
#' @importFrom Matrix Diagonal kronecker t
#' @seealso \code{\link[mgcv]{plot.gam}}, \code{\link[mgcv]{predict.gam}} which this routine is
#'   based on.
#' @author Fabian Scheipl
coef.pffr <- function(object, raw=FALSE, se=TRUE, freq=FALSE, sandwich=FALSE,
                      seWithMean=TRUE, n1=100, n2=40, n3=20, Ktt=NULL, ...){
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
        x <- get.var(trm$term[1], object$model)
        if(trm$dim==1) {
          xg <- if(is.factor(x)) {
            unique(x)
          } else seq(min(x), max(x), length=n1)
          d <- data.frame(xg)
          colnames(d) <- trm$term
          attr(d, "xm") <- xg
        }
        if(is.pcre) {
          ng <- n2
          xg <- if(is.factor(x)) {
            unique(x)
          } else seq(min(x), max(x),length=ng)

          which.pcre <- which(sapply(object$pffr$pcreterms, `[[`, "idname")
                              == trm$term[1])
          pcreterm <- object$pffr$pcreterms[[which.pcre]]

          yg <- seq(min(pcreterm$yind), max(pcreterm$yind), l=ng)

          # interpolate given eigenfunctions to grid values:
          efcts.grid <- sapply(colnames(pcreterm$efunctions),
                               function(nm){
                                 approx(x=pcreterm$yind,
                                        y=pcreterm$efunctions[, nm],
                                        xout=yg,
                                        method = "linear")$y
                               })
          efcts.grid <- data.frame(efcts.grid[rep(1:ng, e=length(xg)),])
          colnames(efcts.grid) <- colnames(pcreterm$efunctions)
          d <- cbind(expand.grid(xg, yg),
                     efcts.grid)
          colnames(d)[1:2] <- c(trm$term[1],
                                paste0(object$pffr$yindname, ".vec"))
          attr(d, "xm") <- xg
          attr(d, "ym") <- yg
        } else {
          if(trm$dim > 1) {
            ng <- ifelse(trm$dim==2, n2, n3)

            varnms <- trm$term

            x <- get.var(trm$term[1], object$model)
            xg <- if(is.factor(x)) {
              unique(x)
            } else seq(min(x), max(x),length=ng)
            y <- get.var(trm$term[2], object$model)
            yg <- if(is.factor(y)) {
              unique(y)
            } else seq(min(y), max(y),length=ng)
            if(length(varnms)==2){
              d <- expand.grid(xg, yg)
              attr(d, "xm") <- xg
              attr(d, "ym") <- yg
            } else {
              z <- get.var(trm$term[3], object$model)
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

        if(is.pcre){
          #sloppy, but effective: temporarily overwrite offending entries
          trm$dim <- 2
          trm$term[2] <- paste0(object$pffr$yindname, ".vec")
        }

        P <- if(trm$dim==1){
          list(x=attr(d, "xm"), xlab=trm$term, xlim=safeRange(attr(d, "xm")))
        } else {
          varnms <-  trm$term
          if(trm$dim==2){
            list(x=attr(d, "xm"), y=attr(d, "ym"), xlab=varnms[1], ylab=varnms[2],
                 ylim=safeRange(attr(d, "ym")), xlim=safeRange(attr(d, "xm")))
          } else {
            if(trm$dim==3){
              list(x=attr(d, "xm"), y=attr(d, "ym"), z=attr(d, "zm"),
                   xlab=varnms[1], ylab=varnms[2], zlab=varnms[3],
                   ylim=safeRange(attr(d, "ym")), xlim=safeRange(attr(d, "xm")),
                   zlim=safeRange(attr(d, "zm")))
            }
          }
        }
        trmind <- trm$first.para:trm$last.para
        P$value <- X%*%object$coefficients[trmind]
        P$coef <- cbind(d, "value"=P$value)
        if(se){
          # use seWithMean if possible:
          if(seWithMean & attr(trm,"nCons")>0){
            cat("using seWithMean for ", trm$label,".\n")
            X1 <- matrix(object$cmX,nrow(X),ncol(object$Vp),byrow=TRUE)
            meanL1 <- trm$meanL1
            if (!is.null(meanL1)) X1 <- X1 / meanL1
            X1[,trmind] <- X
            P$se <- sqrt(rowSums((X1%*%covmat)*X1))
          } else {
            P$se <- sqrt(rowSums((X%*%covmat[trmind, trmind])*X))
          }
          P$coef <- cbind(P$coef, se=P$se)
        }

        P$dim <- trm$dim
        return(P)
      }

      trm <- object$smooth[[i]]
      is.pcre <- "pcre.random.effect" %in% class(trm)

      #FIXME: this fails for pcre-terms with >2 FPCs...!
      if(trm$dim > 3 && !is.pcre){
        warning("can't deal with smooths with more than 3 dimensions, returning NULL for ",
                shrtlbls[names(object$smooth)[i] == unlist(object$pffr$labelmap)])
        return(NULL)
      }

      d <- makeDataGrid(trm)
      P <- getP(trm, d)

      #browser()
      # get proper labeling
      P$main <- shrtlbls[names(object$smooth)[i] == unlist(object$pffr$labelmap)]
      which <- match(names(object$smooth)[i], object$pffr$labelmap)
      if(which %in% object$pffr$where$ff){
        which.ff <- which(object$pffr$where$ff == which)
        P$ylab <- object$pffr$yindname
        xlab <- deparse(as.call(formula(paste("~",names(object$pffr$ff)[which.ff]))[[2]])$xind)
        if(xlab=="NULL") xlab <- "xindex"
        P$xlab <- xlab
      }
      if(which %in% object$pffr$where$sff){
        which.sff <- which(object$pffr$where$sff == which)
        P$ylab <- object$pffr$yindname
        xlab <- deparse(as.call(formula(paste("~",names(object$pffr$ff)[which.sff]))[[2]])$xind)
        if(xlab=="NULL") xlab <- "xindex"
        P$xlab <- xlab
        P$zlab <- gsub(".mat$", "", object$pffr$ff[[which.sff]]$xname)
      }

      return(P)
    }

    bread <- if(freq){
      object$Ve
    } else {
      object$Vp
    }
    if(sandwich){
      X <- predict(object, type = "lpmatrix", reformat=FALSE)
      bread <- bread/object$sig2
      res <- residuals(object)
      if(is.null(Ktt)){
        # get estimate of Cov(eps_i(t)) = K(t,t')
        # stopifnot(require(Matrix))
        Ktt <- Reduce("+",  lapply(1:nrow(res), function(i) tcrossprod(res[i,])))/nrow(res)
      }
      #Chen/Wang, Sec. 2.1: M = X' V^-1 (Y-eta)(Y-eta)' V^-1 X  with V ^-1 = diag(sigma^-2)
      # since the estimate is under working independence....
      meat <- (t(X)%*%kronecker(Diagonal(nrow(res)), Ktt))%*%X / (object$scale^2)
      covmat <- as.matrix(bread %*% meat %*% bread)
    } else {
      covmat <- bread
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

#' Summary for a pffr fit
#'
#' Take a fitted \code{pffr}-object and produce summaries from it.
#' See \code{\link[mgcv]{summary.gam}()} for details.
#'
#' @param object a fitted \code{pffr}-object
#' @param ... see \code{\link[mgcv]{summary.gam}()} for options.
#'
#' @return A list with summary information, see \code{\link[mgcv]{summary.gam}()}
#' @export
#' @method summary pffr
#' @importFrom mgcv summary.gam
#' @author Fabian Scheipl, adapted from \code{\link[mgcv]{summary.gam}()} by Simon Wood, Henric Nilsson
summary.pffr <- function (object, ...) {
  call <- match.call()
  call[[1]] <- mgcv::summary.gam
  ## drop "pffr" class and replace <object> with changed value s.t. method dispatch works without glitches
  ## if we don't do this, summary.gam will call predict on the object if n>3000 & freq==TRUE
  ## and this predict-call gets dispatched to predict.pffr which dispatches back
  ## to predict.gam. Somewhere along the way an index variable get's lost and
  ## shit breaks down.
  class(object) <- class(object)[!(class(object) %in% "pffr")]
  call$object <- as.name("object")
  ret <- eval(call)

  ret$formula <- object$pffr$formula

  # make short labels for display
  shrtlbls <- getShrtlbls(object)

  if(!is.null(ret$s.table)){
    rownames(ret$s.table) <- sapply(rownames(ret$s.table),
                                    function(x){
                                      shrtlbls[pmatch(x, unlist(object$pffr$labelmap))]
                                    })
  }
  class(ret) <- c("summary.pffr", class(ret))
  if(!object$pffr$sparseOrNongrid) {
    ret$n  <- paste(ret$n, " (", object$pffr$nobs," x ", object$pffr$nyindex, ")", sep="")
  } else {
    ret$n  <- paste(ret$n, " (in ", object$pffr$nobs," curves)", sep="")
  }
  return(ret)
}

#' Print method for summary of a pffr fit
#'
#' Pretty printing for a \code{summary.pffr}-object.
#' See \code{\link[mgcv]{print.summary.gam}()} for details.
#'
#' @param x a fitted \code{pffr}-object
#' @param digits controls number of digits printed in output.
#' @param signif.stars Should significance stars be printed alongside output?
#' @param ... not used
#'
#' @return A \code{\link{summary.pffr}} object
#' @method print summary.pffr
#' @importFrom stats printCoefmat
#' @export
#' @author Fabian Scheipl, adapted from \code{\link[mgcv]{print.summary.gam}()} by Simon Wood, Henric Nilsson
print.summary.pffr <- function(x, digits = max(3, getOption("digits") - 3),
                               signif.stars = getOption("show.signif.stars"), ...){
  # mostly identical to print.summary.gam
  print(x$family)
  cat("Formula:\n")
  print(x$formula)
  if (length(x$p.coeff)>0)
  { cat("\nConstant coefficients:\n")
    printCoefmat(x$p.table, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
  }
  cat("\n")
  if(x$m>0)
  { cat("Smooth terms & functional coefficients:\n")
    printCoefmat(x$s.table, digits = digits, signif.stars = signif.stars, has.Pvalue = TRUE, na.print = "NA",cs.ind=1, ...)
  }
  cat("\nR-sq.(adj) = ",formatC(x$r.sq,digits=3,width=5))
  if (length(x$dev.expl)>0) cat("   Deviance explained = ",formatC(x$dev.expl*100,digits=3,width=4),"%\n",sep="")

  if (!is.null(x$method)&&!(x$method%in%c("PQL","lme.ML","lme.REML")))
    cat(x$method," score = ",formatC(x$sp.criterion,digits=5),sep="")

  cat("  Scale est. = ",formatC(x$scale,digits=5,width=8,flag="-"),"  n = ",x$n,"\n",sep="")
  invisible(x)
}

#' QQ plots for pffr model residuals
#'
#' This is simply a wrapper for \code{\link[mgcv]{qq.gam}()}.
#'
#' @param object a fitted \code{\link{pffr}}-object
#' @inheritParams mgcv::qq.gam
#' @importFrom mgcv qq.gam
#' @export
qq.pffr <- function (object, rep = 0, level = 0.9, s.rep = 10, type = c("deviance",
                                                                        "pearson", "response"), pch = ".", rl.col = 2, rep.col = "gray80",
                     ...) {
  if (!inherits(object, "pffr"))
    stop("`object' is not of class \"pffr\"")
  call <- match.call()
  # drop pffr-class so only gam-methods are used on object
  class(object) <- class(object)[-1]
  call$object <- object
  call[[1]] <- mgcv::qq.gam
  eval(call)
}

#' Some diagnostics for a fitted pffr model
#'
#' This is simply a wrapper for \code{\link[mgcv]{gam.check}()}.
#'
#' @param b a fitted \code{\link{pffr}}-object
#' @inheritParams mgcv::gam.check
#' @param rep passed to \code{\link[mgcv]{qq.gam}} when \code{old.style} is \code{FALSE}.
#' @param level passed to \code{\link[mgcv]{qq.gam}} when \code{old.style} is \code{FALSE}.
#' @param rl.col passed to \code{\link[mgcv]{qq.gam}} when \code{old.style} is \code{FALSE}.
#' @param rep.col passed to \code{\link[mgcv]{qq.gam}} when \code{old.style} is \code{FALSE}.
#' @export
pffr.check <- function (b, old.style = FALSE, type = c("deviance", "pearson",
                                                       "response"), k.sample = 5000, k.rep = 200, rep = 0, level = 0.9,
                        rl.col = 2, rep.col = "gray80", ...)  {
  if (!inherits(b, "pffr"))
    stop("`object' is not of class \"pffr\"")
  call <- match.call()
  # drop pffr-class so only gam-methods are used on b
  class(b) <- class(b)[-1]
  call$b <- b
  call[[1]] <- mgcv::gam.check
  eval(call)
}

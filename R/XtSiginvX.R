#' Internal computation function
#'
#' Internal function used compute the products
#' in cross-sectional VB algorithm and Gibbs sampler
#'
#' @param tx transpose of the X design matrix
#' @param siginv inverse variance matrix
#' @param y outcome matrix. if \code{NULL}, function computes
#' first product; if not, function computes second product.
#'
#' @author Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
#'
Xt_siginv_X = function(tx, siginv, y = NULL){

  D = dim(siginv)[1]
  I = dim(tx)[2] / D
  if(is.null(y)){
    ret.mat = matrix(0, nrow = dim(tx)[1], ncol = dim(tx)[1])
    for(i in 1:I){
      ind.cur = (D * (i - 1) + 1) : (D*i)
      prod.cur = tx[,ind.cur] %*% siginv %*% t(tx[,ind.cur])
      ret.mat = ret.mat + prod.cur
    }
  } else if(!is.null(y)){
    ret.mat = matrix(0, nrow = dim(tx)[1], ncol = 1)
    for(i in 1:I){
      ind.cur = (D * (i - 1) + 1) : (D*i)
      prod.cur = tx[,ind.cur] %*% siginv %*% y[ind.cur]
      ret.mat = ret.mat + prod.cur
    }
  }
  return(ret.mat)
}

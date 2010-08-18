\name{fosr}
\Rdversion{1.1}
\alias{fosr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function-on-scalar regression}
\description{
  Fit linear regression with functional responses and scalar predictors,
  with efficient selection of optimal smoothing parameters.
}
\usage{
fosr(fdobj, Z, L = NULL, eval.pts = seq(min(fdobj$basis$range), 
max(fdobj$basis$range), length.out = 201), method = "OLS", 
gam.method = "REML", lambda = NULL, multi.sp = FALSE, max.iter = 1, 
maxlam = NULL, cv1 = FALSE, scale = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fdobj}{
a functional data object (class \code{fd}) giving the functional responses.}
  \item{Z}{
the model matrix, whose columns represent scalar predictors.}
  \item{L}{
a row vector or matrix of linear contrasts of the coefficient functions, to be restricted to equal zero.}
  \item{eval.pts}{
argument values at which the coefficient functions will be evaluated.}
  \item{method}{
estimation method: "OLS" for penalized ordinary least squares, "GLS" for penalized generalized least squares.}
  \item{gam.method}{
smoothing parameter selection method, to be passed to \code{gam()}: "REML" for restricted maximum likelihood, "GCV.Cp" for generalized cross-validation.}
  \item{lambda}{smoothing parameter value.  If \code{NULL}, the smoothing parameter(s) will be estimated.  See Details.}
  \item{multi.sp}{
a logical value indicating whether separate smoothing parameters should be estimated for each coefficient function.  Currently must be \code{FALSE} if \code{method = "OLS"}.}
  \item{max.iter}{
maximum number of iterations if \code{method = "GLS"}.}
  \item{maxlam}{
maximum smoothing parameter value to consider (when \code{lamvec=NULL}; see \code{\link{lofocv}})}
  \item{cv1}{
logical value indicating whether a cross-validation score should be computed even if a fixed\code{lambda} is specified (when \code{method = "OLS"}).}
  \item{scale}{
logical value or vector determining scaling of the matrix \code{Z} (see \code{\link{scale}}, to which the value of this argument is passed).}
}
\details{There are three types of values for argument \code{lambda}:\enumerate{\item if \code{NULL}, the smoothing parameter is estimated by \code{gam} (package \code{mgcv}) if \code{method = "GLS"}, or by \code{optimize} if \code{method = "GLS"};
\item if a scalar, this value is used as the smoothing parameter (but only for the initial model, if \code{method = "GLS"});
\item if a vector, this is used as a grid of values for optimizing the cross-validation score (provided \code{method = "OLS"}; otherwise an error message is issued).}

Please note that currently, if \code{multi.sp = TRUE}, then \code{lambda} must be \code{NULL} and \code{method} must be \code{"GLS"}.
}
\value{
An object of class \code{fosr}, which is a list with the following elements:
\item{B}{matrix of basis coefficients for the estimated coefficient functions.}
\item{yhat}{an object of class \code{fd} giving the fitted values for the functional responses.}
\item{est.func}{matrix of values of the coefficient function estimates at the points given by \code{eval.pts}.}
\item{se.func}{matrix of values of the standard error estimates for the coefficient functions, at the points given by \code{eval.pts}.}
\item{eval.pts}{points at which the coefficient functions are evaluated.}
\item{fit}{fit object outputted by \code{amc()}.}
\item{edf}{effective degrees of freedom of the fit.}
\item{lambda}{smoothing parameter, or vector of smoothing parameters.}
\item{cv}{cross-validated integrated squared error if \code{method="OLS"}, otherwise \code{NULL}.}
\item{roughness}{value of the roughness penalty.}
}
\references{
Ramsay, J. O., and Silverman, B. W. (2005).  \emph{Functional Data Analysis}, 2nd ed., Chapter 13.  New York: Springer.

Reiss, P. T., Huang, L., and Mennes, M. (2010).  Fast function-on-scalar regression with penalized basis expansions.  \emph{International Journal of Biostatistics}, to appear.  Available at \url{http://works.bepress.com/phil_reiss/16/}}
\author{
Philip Reiss \email{phil.reiss@nyumc.org}}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.fosr}}
}
\examples{
# The first two lines, adapted from help(fRegress) in package fda,
# set up a functional data object representing daily average 
# temperatures at 35 sites in Canada
daybasis25 <- create.fourier.basis(rangeval=c(0, 365), nbasis=25,
                  axes=list('axesIntervals'))
Temp.fd <- with(CanadianWeather, smooth.basisPar(day.5,
                dailyAv[,,'Temperature.C'], daybasis25)$fd)

modmat = cbind(1, model.matrix(~ factor(CanadianWeather$region) - 1))
constraints = matrix(c(0,1,1,1,1), 1)

# Penalized OLS with smoothing parameter chosen by grid search
olsmod = fosr(Temp.fd, modmat, constraints, method="OLS", lambda=100*10:30)
plot(olsmod, 1)

# Penalized GLS
glsmod = fosr(Temp.fd, modmat, constraints, method="GLS")
plot(glsmod, 1)
}

\keyword{models}
\keyword{regression}
\keyword{smooth}

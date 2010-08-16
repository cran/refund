\name{fpcr}
\Rdversion{1.1}
\alias{fpcr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functional principal component regression}
\description{
  Implements functional principal component regression (Reiss and Ogden, 2007, 2010) for generalized linear models with scalar responses and functional predictors.}
\usage{
fpcr(y, sigmat = NULL, basismat = NULL, halfpens = NULL, fdobj = NULL, 
     eval.pts = if (!is.null(fdobj)) seq(fdobj$basis$rangeval[1],
                    fdobj$basis$rangeval[2], length.out=401) else NULL,
     nc, covt = NULL, mean.signal.term = FALSE, 
     family = gaussian(), method = "REML", gamma = 1, sp = NULL, 
     min.sp = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
scalar outcome vector.}
  \item{sigmat}{
matrix of signals: n x N, where n is the length of y and N is the
 number of sites at which each signal is defined.}
  \item{basismat}{
N x K matrix of values of K basis functions at the N sites}
  \item{halfpens}{
a list, each component of which is a matrix M with K columns such 
that M'M defines a penalty on the basis coefficients.}
\item{fdobj}{functional data object (class \code{fd}) giving the functional predictors.}
\item{eval.pts}{points at which to evaluate the signals and the coefficient function.}
  \item{nc}{
number of principal components.}
  \item{covt}{
covariates: an n-row matrix, or a vector of length n.}
  \item{mean.signal.term}{
logical: should the mean of each subject's signal be 
included as a covariate?}
  \item{family}{
generalized linear model family.}
  \item{method, gamma, sp, min.sp}{
arguments passed to function \code{gam} in package \code{mgcv}; see the \code{gam} documentation for details.}
}
\details{
The functional predictors and penalty matrix are given either explicitly using \code{sigmat}, \code{basismat} and \code{halfpens}, or in functional data form (see the \code{fda} package of Ramsay, Hooker and Graves, 2009) using arguments \code{fdobj} and \code{eval.pts}.  If the latter two arguments are non-null, the values of the former three are ignored.  Optimal smoothing is achieved by means of the \code{gam} function in package \code{mgcv} (Wood, 2006).}
\value{
An object of class \code{gam} (see \code{gamObject} in the \code{mgcv} package documentation), with two 
additional components: \code{nc} (same as argument \code{nc}) and \code{fhat} (coefficient 
 function estimate).}
\references{
Ramsay, J. O., Hooker, G., and Graves, S. (2009).
\emph{Functional Data Analysis with R and MATLAB}.
New York: Springer.

Reiss, P. T., and Ogden, R. T. (2007).  Functional principal component regression and functional partial least squares.  \emph{Journal of the American Statistical Association}, 102, 984--996.

Reiss, P. T., and Ogden, R. T. (2010).  Functional generalized linear models with images as predictors.  \emph{Biometrics}, 66, 61--69.

Wood, S. N. (2006). \emph{Generalized Additive Models: An Introduction with R}. Boca Raton, FL: Chapman &
  Hall.}

\author{
Philip Reiss \email{phil.reiss@nyumc.org}}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##### Octane data example #####
require(pls)
data(gasoline)

# Create the requisite functional data objects
bbasis = create.bspline.basis(c(900, 1700), 40)
wavelengths = 2*450:850
nir = matrix(NA, 401, 60)
for (i in 1:60) nir[ , i] = gasoline$NIR[i, ]
# Why not just take transpose of gasoline$NIR above?  
# Because for some reason it leads to an error in the next statement
gas.fd = smooth.basisPar(wavelengths, nir, bbasis)$fd

# Method 1: Call fpcr with fdobj argument
gasmod1 = fpcr(gasoline$octane, fdobj = gas.fd, nc = 30)

# Method 2: Call fpcr with explicit signal, basis, and penalty matrices
# --the following lines use functions from package 'fda' to derive
# these matrices so as to yield identical results to the above
mybasismat = eval.basis(wavelengths, bbasis)
bb2 = quadset(basisobj = bbasis)
myhalfpen = crossprod(diag(sqrt(bb2$quadvals[ , 2])), eval.basis(bb2, bb2$quadvals[ , 1], 2))

gasmod2 = fpcr(gasoline$octane, sigmat = gasoline$NIR, basismat = mybasismat, 
               halfpens = list(myhalfpen), nc = 30)

# Check that both calls yield identical estimates
plot(wavelengths, gasmod1$fhat, type='l', xlab='Wavelength', 
     ylab='Coefficient function')
lines(wavelengths, gasmod2$fhat, col=2)
}

\keyword{models}
\keyword{regression}
\keyword{smooth}

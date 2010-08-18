\name{amc}
\Rdversion{1.1}
\alias{amc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Additive model with constraints}
\description{
An internal function, called by \code{fosr()}, that fits additive models with linear constraints via a call to the \code{gam()}
function in the \code{mgcv} package.}
\usage{
amc(y, Xmat, S, gam.method = "REML", C = NULL, lambda = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
response vector.}
  \item{Xmat}{
design matrix.}
  \item{S}{
list of penalty matrices.}
  \item{gam.method}{
smoothing parameter selection method: "REML" for restricted maximum likelihood, "GCV.Cp" for generalized cross-validation.}
  \item{C}{
matrix of linear constraints.  Dimension should be number of constraints times \code{ncol(Xmat)}.}
  \item{lambda}{
smoothing parameter value.  If \code{NULL}, the smoothing parameter(s) will be estimated.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with the following elements:
\item{gam}{the \code{gam} object returned by \code{gam()}.}
\item{coefficients}{coefficients with respect to design matrix \code{Xmat}, derived from the \code{gam()} fit.}
\item{Vp, GinvXt}{outputs used by \code{fosr()}.}
\item{method}{the \code{gam.method} argument of the call to \code{amc()}.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Philip Reiss \email{phil.reiss@nyumc.org}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fosr}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
\keyword{smooth}

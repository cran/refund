% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fpc.R
\name{smooth.construct.fpc.smooth.spec}
\alias{smooth.construct.fpc.smooth.spec}
\title{Basis constructor for FPC terms}
\usage{
\method{smooth.construct}{fpc.smooth.spec}(object, data, knots)
}
\arguments{
\item{object}{a \code{fpc.smooth.spec} object, usually generated by a
term \code{s(x, bs="fpc")}; see Details.}

\item{data}{a list containing the data (including any \code{by} variable)
required by this term, with names corresponding to \code{object$term}
(and \code{object$by}). Only the first element of this list is used.}

\item{knots}{not used, but required by the generic \code{smooth.construct}.}
}
\value{
An object of class \code{"fpc.smooth"}. In addtional to the elements
  listed in \code{{smooth.construct}}, the object will contain
  \item{sm}{the smooth that is fit in order to generate the basis matrix
    over \code{object$term}}
  \item{V.A}{the matrix of principal components}
}
\description{
Basis constructor for FPC terms
}
\details{
\code{object} must contain an \code{xt} element. This is a list that can
  contain the following elements:
\describe{
  \item{X}{(required) matrix of functional predictors}
  \item{method}{(required) the method of finding principal components;
    options include \code{"svd"} (unconstrained), \code{"fpca.sc"},
    \code{"fpca.face"}, or \code{"fpca.ssvd"}}
  \item{npc}{(optional) the number of PC's to retain}
  \item{pve}{(only needed if \code{npc} not supplied) the percent variance
    explained used to determine \code{npc}}
  \item{penalize}{(required) if \code{FALSE}, the smoothing parameter is
    set to 0}
  \item{bs}{the basis class used to pre-smooth \code{X}; default is \code{"ps"}}
}

Any additional options for the pre-smoothing basis (e.g. \code{k}, \code{m},
  etc.) can be supplied in the corresponding elements of \code{object}.
  See \code{[mgcv]{s}} for a full list of options.
}
\references{
Reiss, P. T., and Ogden, R. T. (2007). Functional principal component
regression and functional partial least squares. \emph{Journal of the
American Statistical Association}, 102, 984-996.
}
\seealso{
\code{{fpcr}}
}
\author{
Jonathan Gellar \email{JGellar@mathematica-mpr.com}
}

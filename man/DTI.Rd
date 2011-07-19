\name{DTI}
\alias{DTI}
\docType{data}
\title{
Diffusion Tensor Imaging: tract profiles and outcomes
}
\description{
Fractional anisotropy tract profiles for the corpus callosum (cca) and the right corticospinal tract (rcst). Accompanying the tract profiles are the subject ID numbers, visit number, total number of scans, multiple sclerosis case status and Paced Auditory Serial Addition Test (pasat) score. We thank Dr. Daniel Reich for making this dataset available.
}
\usage{data(DTI)}
\format{
  A data frame made up of
  \describe{
    \item{\code{cca} }{A 382 x 93 matrix of fractional anisotropy tract profiles from the corpus callosum;}
    \item{\code{rcst} }{A 382 x 55 matrix of fractional anisotropy tract profiles from the right corticospinal tract;}
    \item{\code{ID} }{Numeric vector of subject ID numbers;}
    \item{\code{visit} }{Numeric vector of the subject-specific visit numbers;}
    \item{\code{Nscans} }{Numeric vector indicating the total number of visits for each subject;}
    \item{\code{case} }{Numeric vector of multiple sclerosis case status: 0 - healthy control, 1 - MS case;}
    \item{\code{pasat }}{Numeric vector containing the PASAT score at each visit.}
  }
}
\references{
Goldsmith, J., Bobb, J., Crainiceanu, C., Caffo, B., and Reich, D. (2011). Penalized Functional Regression. \emph{Journal of Computational and Graphical Statistics}, .., .. - ..

Goldsmith, J., Crainiceanu, C., Caffo, B., and Reich, D. (2010). Longitudinal Penalized Functional Regression. www.bepress.com/jhubiostat/paper216
}


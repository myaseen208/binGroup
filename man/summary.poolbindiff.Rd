\name{summary.poolbindiff}
\alias{summary.poolbin}
\alias{summary.poolbindiff}
\title{Summary methods for "poolbin" and "poolbindiff"}
\description{
Summary Method for One-Sample and Two-sample confidence intervals (various pool sizes)
}
\usage{
\method{summary}{poolbin}(object, scale = object$scale, ...)
\method{summary}{poolbindiff}(object, scale = object$scale, ...)
}

\arguments{
  \item{object}{ An object of class "poolbin" or "poolbindiff" (\code{\link{pooledBin}}, \code{\link{pooledBinDiff}})}
  \item{scale}{ A coefficient to scale the point estimate and interval bounds}
  \item{\dots}{ further arguments to be passed to \code{print}}
}

\value{
A print out
}

\keyword{ print }
\name{print.poolbindiff}
\alias{print.poolbin}
\alias{print.poolbindiff}

\title{
Print methods for classes "poolbin" and "poolbindiff"
}
\description{
Print methods for objects of classes \code{"poolbin"} and \code{"poolbindiff"}
}
\usage{
\method{print}{poolbin}(x, scale = x$scale, ...)
\method{print}{poolbindiff}(x, scale = x$scale, ...)
}
\arguments{
  \item{x}{ An object of class "poolbin" or "poolbindiff" (\code{\link{pooledBin}}, \code{\link{pooledBinDiff}})}
  \item{scale}{ A coefficient to scale the point estimate and interval bounds }
  \item{\dots}{ further arguments to be passed to \code{print}}
}

\value{
A print out.
}

\keyword{ print }

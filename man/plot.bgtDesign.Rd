\name{plot.bgtDesign}
\alias{plot.nDesign}
\alias{plot.sDesign}

\title{Plot Results of nDesign or sDesign}

\description{
Plots the results of the iteration performed by \kbd{nDesign} or \kbd{sDesign} in order to find an experimental design with sufficient power.
}
\usage{
\method{plot}{nDesign}(x, \dots)
\method{plot}{sDesign}(x, \dots)
}

\arguments{
  \item{x}{ An object of class \kbd{"nDesign"} or \kbd{"sDesign"} as can be created by \code{\link{nDesign}} and \code{\link{sDesign}}, respectively.}
  \item{\dots}{ further arguments to be passed to \code{\link{plot}} }
}

\examples{

plot(nDesign(nmax=100, s=30, delta=0.004, p.hyp=0.005,
 alternative="less", method="Score", power=0.8))

}
\keyword{ hplot }
\concept{ group testing}
\concept{ sample size }
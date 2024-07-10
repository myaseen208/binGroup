\name{plot.poolbin}
\alias{plot.poolbin}

\title{
Diagnostic line fit for pool.bin objects
}
\description{
Diagnostic line fit for pool.bin objects
}
\usage{
\method{plot}{poolbin}(x, pch = 16, refline = TRUE, printR2 = TRUE, ...)
}

\arguments{
  \item{x}{
an object of class \code{"poolbin"}, as can be created by \code{\link{pooledBin}}
}
  \item{pch}{symbol type for plotting, see \code{?pch}
}
  \item{refline}{
logical, whether or not to add a reference line to the plot
}
  \item{printR2}{
logical, whether to print R^2 
}
  \item{\dots}{
further arguments passed to \code{plot}
}
}
\details{
Brad: Maybe drop a line here...
}
\value{
a plot
}
\references{
Brad: Any reference for this?
}
\author{
Brad Biggerstaff
}

\examples{

test <- pooledBin(x=c(8,3,1,0), m=c(50,20,10,5),
n=c(10,10,10,10), pt.method="firth", ci.method="lrt")

plot(test)
}

\keyword{ hplot }

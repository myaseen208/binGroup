\name{print.summary.gt}
\Rdversion{1.1}
\alias{print.summary.gt.mp}
\alias{print.summary.gt}

\title{
Print Functions for summary.gt.mp and summary.gt
}
\description{
Print function for objects obtained by calling \code{\link{summary.gt.mp}} and \code{\link{summary.gt}}
}
\usage{
\method{print}{summary.gt.mp}(x, digits = max(3, getOption("digits") - 3),
 signif.stars = getOption("show.signif.stars"), ...)
\method{print}{summary.gt}(x, digits = max(3, getOption("digits") - 3),
 signif.stars = getOption("show.signif.stars"), ...)
}

\arguments{
  \item{x}{
an object of class \kbd{"summary.gt.mp"} or \kbd{"summary.gt"}
}
  \item{digits}{
digits for rounding
}
  \item{signif.stars}{
logical, indicating whether significance stars should be shown
}
  \item{\dots}{
 further arguments to be passed to \kbd{printCoefmat}
}
}

\value{
A print out.
}

\keyword{print}
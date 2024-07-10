\name{print.bgt}
\alias{print.bgtCI}
\alias{print.binCI}
\alias{print.bgtTest}
\alias{print.binTest}
\alias{print.bgtvs}

\title{ Print Functions for Group Testing CIs and Tests for One Proportion }

\description{ Print objects of class \kbd{"bgtCI"}, \kbd{"bgtTest"}, \kbd{"bgtvs"}, \kbd{"binCI"}, and \kbd{"binTest"}  }
\usage{
\method{print}{bgtCI}(x, \dots)
\method{print}{binCI}(x, \dots)
\method{print}{bgtTest}(x, \dots)
\method{print}{binTest}(x, \dots)
\method{print}{bgtvs}(x, \dots)
}

\arguments{
  \item{x}{ an object of the corresponding class}
  \item{\dots}{currently only \code{digits} is passed to \code{signif} for appropriate rounding}
}

\value{
A print out.
}

\keyword{ print }
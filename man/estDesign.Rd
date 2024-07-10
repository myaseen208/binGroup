\name{estDesign}
\alias{estDesign}
\alias{msep}
\title{Sample Size Iteration Depending on Minimal MSE in One-Parameter Group Testing}
\description{
Find the group size s for a fixed number of assays n and an assumed true proportion p.tr for which
the mean squared error (mse) of the point estimator is minimal and bias is within a restriction. For experimental design in binomial group testing
as recommended by Swallow (1985), if main objective is estimation.
}
\usage{
estDesign(n, smax, p.tr, biasrest = 0.05)
}
\arguments{
  \item{n}{ integer, fixed sample size (number of assays) }
  \item{smax}{ integer, maximal group size allowed in planning of the design }
  \item{p.tr}{ assumed true proportion of the 'positive' trait in the population to be tested, specify as a value between 0 and 1 }
  \item{biasrest}{ value between 0 and 1 specifying the absolute bias maximally allowed }
}

\details{ Swallow (1985) recommends to use the upper bound of the expected range of true proportion p.tr for optimization of tzhe design.
For further details see the reference.
Up to now, specify n<1020.
}
\value{the group size s, for which the mse of the estimator is minimal for the given n, p.tr
 or the group size s for which bias restriction biasrest is just not violated,
 and for this particular group size s: a list containing:
  \item{varp }{the variance of the estimator }
  \item{mse }{the mean square error of the estimator } 
  \item{bias }{the bias of the estimator }
  \item{exp }{the expected value of the estimator }

}
\references{ \emph{Swallow WH (1985)} Group testing for estimating infection rates and probabilities of disease transmission. Phytopathology Vol.75, N.8, 882-889.
  }

\author{
Frank Schaarschmidt
}

\seealso{ \code{\link{nDesign}}, \code{\link{sDesign}} for choice of the binomial group testing design according to the power in a hypothesis test }
\examples{

### Compare table 1 in Swallow(1985),885:

estDesign(n=10, smax=100, p.tr=0.001)

estDesign(n=10, smax=100, p.tr=0.01)

estDesign(n=25, smax=100, p.tr=0.05)

estDesign(n=40, smax=100, p.tr=0.25)

estDesign(n=200, smax=100, p.tr=0.3)


}
\keyword{ htest }
\concept{ group testing}
\concept{ sample size}



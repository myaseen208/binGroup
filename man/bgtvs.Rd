\name{bgtvs}
\alias{bgtvs}
\title{ Confidence Interval for One Proportion in Group Testing with Variable Group Sizes }
\description{
Calculates confidence intervals for a single proportion in binomial group testing
if groups of different size are evaluated
}
\usage{
bgtvs(n, s, y, conf.level = 0.95, alternative = "two.sided",
 maxiter = 100)
}
\arguments{
  \item{n}{vector of integer values, specifying the number of groups of the different sizes }
  \item{s}{vector of integer values, specifying the group sizes, must be of same length as \code{n} }
  \item{y}{vector of integer values, specifying the number of positive groups among the \code{n} groups tested}
  \item{conf.level}{a single numeric value, the confidence level of the interval }
  \item{alternative}{a character string, with options "two.sided", "less", "greater" }
  \item{maxiter}{maximal number steps in iteration of confidence limits }
}
\details{
 Hepworth (1996) describes methods for constructing confidence intervals in binomial group testing,
if groups of different size are used. Currently, only the exact method (Hepworth, 1996, equation5, Table.5)
is implemented. Note, that the algorithm becomes computationally very expensive if the number of different
groups becomes larger than 3.

}
\value{A list containing
  \item{conf.int}{a numeric vector, the lower and upper limits of the confidence interval}
  \item{estimate}{the point estimate}
  \item{conf.level}{as input}
  \item{alternative}{as imput}
  \item{input}{a matrix containing the input values of n (number of groups), s (group size), and y (number of positive pools)}

  moreover, some of the input arguments.
}
\references{ \emph{Hepworth, G (1996):} Exact confidence intervals for proportions estimated by group testing. \emph{Biometrics 52, 1134-1146.}
 }


\author{
Frank Schaarschmidt
}

\seealso{ \code{\link{pooledBin}} for asymptotic methods to calculate confidence intervals for one proportion in designs with a number of different pool sizes. Note that \code{pooledBin} can handle larger number of different pool sizes than \code{bgtvs} }


\examples{
# Consider a very simple example,
# given in Hepworth (1996), table 5:
# 2 groups each containing 5 units,
# and 3 groups, each containing 2 units

# In the first setting (n=2, s=5) y=1 positive group 
# has been observed, in the second setting (n=3, s=2),
# y=2 positive have been observed.

bgtvs(n=c(2,3), s=c(5,2), y=c(1,2)) 

###############################################

# Recalculate the example given in
# Hepworth (1996), table 5:

 bgtvs(n=c(2,3), s=c(5,2), y=c(0,0))
 bgtvs(n=c(2,3), s=c(5,2), y=c(0,1))
 bgtvs(n=c(2,3), s=c(5,2), y=c(0,2))
 bgtvs(n=c(2,3), s=c(5,2), y=c(0,3))
 bgtvs(n=c(2,3), s=c(5,2), y=c(1,0))
 bgtvs(n=c(2,3), s=c(5,2), y=c(1,1))
 bgtvs(n=c(2,3), s=c(5,2), y=c(1,2))
 bgtvs(n=c(2,3), s=c(5,2), y=c(1,3))
 bgtvs(n=c(2,3), s=c(5,2), y=c(2,0))
 bgtvs(n=c(2,3), s=c(5,2), y=c(2,1))
 bgtvs(n=c(2,3), s=c(5,2), y=c(2,2))
 bgtvs(n=c(2,3), s=c(5,2), y=c(2,3))

}

\keyword{ htest }
\concept{group testing}
\concept{variable group size}

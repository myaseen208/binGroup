\name{binCI}
\alias{binCI}
\alias{binCP}
\alias{binBlaker}
\alias{binSOC}
\alias{binAC}
\alias{binWilson}
\alias{binWald}

\title{ Confidence Intervals for One Binomial Proportion }
\description{
Calculates the exact Clopper-Pearson and Blaker, the asymptotic second-order corrected,
Wilson, Agresti-Coull and Wald confidence interval for a single binomial proportion
}
\usage{
binCI(n, y, conf.level = 0.95, alternative = "two.sided",
 method = "CP")

binCP(n, y, conf.level=0.95, alternative="two.sided")
binBlaker(n,y,conf.level=0.95, tolerance=1e-04, alternative="two.sided")
binAC(n, y, conf.level=0.95, alternative="two.sided")
binSOC(n, y,conf.level=0.95,alternative="two.sided")
binWald(n, y, conf.level=0.95, alternative="two.sided")
binWilson(n, y,conf.level=0.95,alternative="two.sided")
}
\arguments{
  \item{n}{ number of trials (number of individuals under observation) }
  \item{y}{ number of successes (number of individuals showing the trait of interest) }
  \item{conf.level}{ nominal confidence level }
  \item{alternative}{ character string defining the alternative hypothesis, either 'two.sided', 'less' or 'greater'
  where 'less' gives the only an upper bound with confidence level=conf.level
  'greater' gives the only a lower bound with confidence level=conf.level
   and 'two.sided' gives a two-sided confidence interval with confidence level=conf.level }
  \item{method}{character string defining the method for CI calculation:
		where "CP" is Clopper-Pearson, an exact tail interval showing symmetric coverage probability (inversion of two one-sided tests),
		"Blaker" is the Blaker interval, an exact interval, inversion of one two.sided test, therefore defined only two.sided, but shorter
                     than the two-sided Clopper-Pearson CI. Both guarantee to contain the true parameter with at least conf.level*100 percent
                     probability,
		"AC" is Agresti-Coull, generalized Agresti-Coull interval, asymptotic method,
		"Score" is Wilson Score, asymptotic method derived from inversion of the Score test,
 		"SOC" is the second order corrected interval, asymptotic method for one-sided problems (for details see Cai, 2005),
		and "Wald" the Wald interval, which cannot be recommended.}
  \item{tolerance}{precision of computation for the bounds of the Blaker interval}
}

\details{
This function allows computation of confidence intervals for a binomial proportion from a standard binomial experiment.
If an actual confidence level greater or equal to that specified in the conf.level argument shall always be guaranteed, the exact method of 
Clopper-Pearson (method="CP") can be recommended for one-sided and the improved method of Blaker (method="Blaker") can be recommended for two-sided
hypotheses. If a mean confidence level close to that specified in the argument conf.level is required, but moderate violation of this level is acceptable,
the Second-Order corrected (method="SOC"), Wilson Score (method="Wilson") or Agresti-Coull (method="AC") might be used, where SOC has the most symmetric
coverage and Wilson and Agresti-Coull are in tendency conservative for the upper bound and proportions close to 0 and for the lower bound and proportions
close to 1. The Wald CI might be used for large number of observations n>10000 or intermediate proportions.

For discussion of CI for a single binomial proportion see Brown et al. (2001) for two-sided and Cai (2005) for one-sided intervals.

}
\value{ A list containing:
  \item{conf.int }{the estimated confidence interval}
  \item{estimate }{the point estimator}

  And the method, conf.level and alternative specified in the function call.
}
\references{ 

\emph{Blaker H (2000)} Confidence curves and improved exact confidence intervals for discrete distributions. \emph{The Canadian Journal of Statistics 28 (4), 783-798.}

\emph{Brown LD, Cai TT, DasGupta A (2001)} Interval estimation for a binomial proportion. \emph{Statistical Science 16 (2), 101-133.}

\emph{Cai TT(2005)} One-sided confidence intervals in discrete distributions. \emph{Journal of Statistical Planning and Inference 131, 63-88.}

}

\author{
Frank Schaarschmidt
}

\seealso{ \code{\link{binom.test}} for the exact confidence interval and test,
\code{\link{binTest}} to calculate p.values of the exact, Score and Wald test.
 }
\examples{

# Default method is the two-sided 95% Clopper-Pearson CI:

binCI(n=200, y=10)

# other methods might result in 
# shorter intervals (but asymetric coverage):

binCI(n=200,y=10, method="Blaker")
binCI(n=200,y=10, method="Score")

}
\keyword{ htest }
\concept{ binomial }
\concept{ confidence interval }
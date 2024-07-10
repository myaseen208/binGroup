\name{bgtCI}

\alias{bgtCI}
\alias{bgtCP}
\alias{bgtBlaker}
\alias{bgtAC}
\alias{bgtWilson}
\alias{bgtWald}
\alias{bgtSOC}

\title{ Confidence Intervals for One Proportion in Binomial Group Testing  }

\description{ Calculates the point estimate, the exact Clopper-Pearson and Blaker CI, the Score test derived 
		Wilson and Agresti-Coull CI, the asymptotic second-order corrected interval fo Cai and the Wald CI for
		a single binomial proportion estimated from a binomial group testing trial. Assumes equal group sizes, 
		an assay method classifying a group as positive if at least one unit in the group is positive, individuals
		units randomly assigned to the groups.
 }
\usage{
bgtCI(n, s, y, conf.level = 0.95,
 alternative = "two.sided", method = "CP")
}

\arguments{
  \item{n}{integer, specifying the number of groups (i.e. assays i.e. observations) }
  \item{s}{integer, specifying the common size of groups i.e. the number of individual units in each group }
  \item{y}{integer, specifying the number of positive groups }
  \item{conf.level}{nominal confidence level of the interval }
  \item{alternative}{character string defining the alternative hypothesis, either 'two.sided', 'less' or 'greater'
  where 'less' gives the only an upper bound with confidence level=conf.level
  'greater' gives the only a lower bound with confidence level=conf.level
   and 'two.sided' gives a two-sided confidence interval with confidence level=conf.level}
   \item{method}{ character string defining the method for CI calculation, where: 
     \kbd{"CP"} is Clopper-Pearson, an exact tail interval showing symmetric coverage probability (inversion of two one-sided tests),
     \kbd{"Blaker"} is the Blaker interval, an exact interval, inversion of one two.sided test, therefore defined only two.sided, but shorter
             than the two-sided Clopper-Pearson CI. Both guarantee to contain the true parameter with at least conf.level*100 percent
             probability,
     \kbd{"AC"} is the Agresti-Coull (generalized Agresti-Coull) interval, asymptotic method,
     \kbd{"Score"} is Wilson Score, asymptotic method derived from inversion of the Score test,
     \kbd{"SOC"} is the second order corrected interval, asymptotic method for one-sided problems (for details see Cai, 2005),
		and \kbd{"Wald"} the Wald interval, which cannot be recommended.
 }

}
\details{

This function allows the computation of confidence intervals for binomial group testing as described in Tebbs & Bilder (2004) and Schaarschmidt (2007).
If an actual confidence level greater or equal to that specified in the conf.level argument shall always be guaranteed, the exact method of 
Clopper-Pearson (\kbd{method="CP"}) can be recommended for one-sided and the improved method of Blaker (2000) (\kbd{method="Blaker"}) can be recommended for two-sided
hypotheses. If a mean confidence level close to that specified in the argument conf.level is required, but moderate violation of this level is acceptable,
the Second-Order corrected (\kbd{method="SOC"}), Wilson Score (\kbd{method="Score"}) or Agresti-Coull (\kbd{method="AC"}) might be used (Brown, Cai, DasGupta, 2001; Cai 2005).}

\value{
A list containing:
  \item{conf.int}{a confidence interval for the proportion}
  \item{estimate}{the point estimator of the proportion}
}

\references{

\emph{Blaker H (2000)}. Confidence curves and improved exact confidence intervals for discrete distributions. \emph{The Canadian Journal of Statistics 28 (4), 783-798.}

\emph{Brown LD, Cai TT, DasGupta A (2001)}. Interval estimation for a binomial proportion. \emph{Statistical Science 16 (2), 101-133.}

\emph{Cai TT (2005)}. One-sided confidence intervals in discrete distributions. \emph{Journal of Statistical Planning and Inference 131, 63-88.}

\emph{Schaarschmidt F (2007).} Experimental design for one-sided confidence intervals or hypothesis tests in binomial group testing. \emph{Communications in Biometry and Crop Science 2 (1), 32-40.} http://agrobiol.sggw.waw.pl/cbcs/

\emph{Tebbs JM & Bilder CR (2004)}. Confidence interval procedures for the probability of disease transmission in multiple-vector-transfer designs. \emph{Journal of Agricultural, Biological and Environmental Statistics, 9 (1), 75-90.}

 }

\author{
Frank Schaarschmidt
}

\seealso{ 
\code{\link{pooledBin}} for asymptotic confidence intervals and \code{\link{bgtvs}} for an exact confidence interval when designs with different group sizes are used

\code{\link{bgtTest}}: for hypothesis tests in binomial group testing
}
\examples{

# See the example in Tebbs and Bilder (2004)
# the two.sided 95-percent 
# Clopper-Pearson as default method:

bgtCI(n=24,y=3,s=7)
bgtCI(n=24,y=3,s=7,conf.level=0.95,
 alternative="two.sided", method="CP")

# other methods:
# Blaker CI is exact but shorter
# than Clopper-Pearson, only two.sided

bgtCI(n=24,y=3,s=7, alternative="two.sided",
 method="Blaker")

# the asymptotic Wilson CI might even
# be shorter: 

bgtCI(n=24,y=3,s=7, alternative="two.sided",
 method="Score")

# one-sided confidence intervals:

bgtCI(n=24,y=3,s=7, alternative="less", method="CP")

# Wilson Score interval is less conservative 
bgtCI(n=24,y=3,s=7, alternative="less", method="Score")

# the second-order corrected CI is even shorter
# in this situation:
bgtCI(n=24,y=3,s=7, alternative="less", method="SOC")

}

\keyword{htest}
\concept{binomial group testing}
\concept{confidence interval}

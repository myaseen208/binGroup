\name{binWidth}
\alias{binWidth}

\title{Expected Confidence Interval Width for One Binomial Proportion }

\description{
Calculation of expected value of the width of confidence intervals in a binomial experiment,
in dependence of the number of trials (number of individuals under observation),  confidence level and an assumed true proportion.
Available for the confidence interval methods in binCI(binGroup).

}
\usage{
binWidth(n, p, conf.level = 0.95,
 alternative = "two.sided", method = "CP")
}
 
\arguments{

  \item{n}{ integer, giving the number of trials (i.e. number of individuals under observation) }
  \item{p}{ assumed true proportion of individuals showing the trait to be estimated }
  \item{conf.level}{ required confidence level of the interval }
  \item{alternative}{ character string, defining the alternative hypothesis, either 'two.sided', 'less' or 'greater'
where 'less' calculates the expected width between the assumed true proportion p and the upper conf.level*100 percent-bound of a one-sided CI,
'greater' calculates the expected width between the assumed true proportion p and the lower conf.level*100 percent-bound of a one-sided CI,
'two.sided' calculates the expected width between the lower and the upper bound of a two-sided conf.level*100 percent-CI.}

  \item{method}{ character string defining the method for CI calculation:
		where \kbd{"CP"} is Clopper-Pearson, an exact tail interval showing symmetric coverage probability (inversion of two one-sided tests),
		\kbd{"Blaker"} is the Blaker interval, an exact interval, inversion of one two.sided test, therefore defined only two.sided, but shorter
                     than the two-sided Clopper-Pearson CI. Both guarantee to contain the true parameter with at least conf.level*100 percent
                     probability,
		\kbd{"AC"} is Agresti-Coull, generalized Agresti-Coull interval, asymptotic method,
		\kbd{"Score"} is Wilson Score, asymptotic method derived from inversion of the Score test,
 		\kbd{"SOC"} is the second order corrected interval, asymptotic method for one-sided problems (for details see Cai, 2005),
		and \kbd{"Wald"} the simple Wald-type interval.
		}

}
\details{
For calculation of expected interval width in the standard binomial estimation see Brown et al. (2001).
}
\value{ A list containing:
  \item{expCIWidth}{the expected value of the width of the confidence interval for the specified arguments}
  
and the alternative, p and n which are specified in the function call.
}


\author{
Frank Schaarschmidt
}

\seealso{ \code{\link{binDesign}} for experimental design for hypothesis testing }
\examples{


# methods differ slightly in length when sample sizes are large:

binWidth(n=200,p=0.02,alternative="two.sided",
 method="CP")$expCIWidth

binWidth(n=200,p=0.02,alternative="two.sided",
 method="Blaker")$expCIWidth

binWidth(n=200,p=0.02,alternative="two.sided",
 method="Score")$expCIWidth

# but do more for small sample sizes and intermediate p:

binWidth(n=20,p=0.2,alternative="two.sided",
 method="CP")$expCIWidth

binWidth(n=20,p=0.2,alternative="two.sided",
 method="Blaker")$expCIWidth

binWidth(n=20,p=0.2,alternative="two.sided",
 method="Score")$expCIWidth

}
\keyword{ htest }
\concept{ sample size }
\concept{ single binomial proportion }
\concept{ interval width }

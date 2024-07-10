\name{binPower}
\alias{binPower}
\alias{binPowerI}

\title{Power Calculation for One Parameter Binomial Problem}

\description{
Closed calculation of the power to reject a hypothesis using confidence intervals for a single 
binomial proportion, for specified sample size n, conf.level and an assumed absolute difference to the
threshold parameter under the null hypothesis. 
}
\usage{
binPower(n, delta, p.hyp, conf.level = 0.95, 
alternative = "two.sided", method = "CP")
}
\arguments{
  \item{n}{number of trials n in the binomial experiment, specify as a single integer  }
  \item{delta}{assumed absolute difference of the true proportion to the threshold proportion under the null hypothesis  }
  \item{p.hyp}{threshold proportion under the null hypothesis }
  \item{conf.level}{nominal confidence level of the interval }
   \item{alternative}{character string, defining the alternative hypothesis, either \kbd{'two.sided'}, \kbd{'less'} or \kbd{'greater'}
where \kbd{'less'} calculates the 'power of the upper confidence limit' for a true proportion \kbd{p.hyp} - \kbd{delta},
\kbd{'greater'} calculates the 'power of the lower confidece limit' for a true proportion of \kbd{p.hyp} + \kbd{delta}.
\kbd{'two.sided'} calculates min(power(\kbd{p.hyp} - \kbd{delta}, \kbd{p.hyp} + \kbd{delta})) for a two-sided CI, thus can result in much lower power.
 }

  \item{method}{ Character string, specifying the confidence interval method (see \code{\link{binCI}}) to be used }

}
\details{
The power of a confidence interval here is defined as the probability that a confidence interval or limit excludes the threshold parameter (\kbd{p.hyp}) of the null hypothesis. 
E.g., the null hypothesis H0: p>= 0.005 can be rejected, if an upper confidence limit for p does not contain p.hyp=0.005.
In case that a delta of 0.002 shall be detectable, this function calculates the probability, that an interval of a given method will
exclude p.hyp=0.005 if the true proportion = 0.003.
Due to discreteness, the power does not increase monotone for increasing sample size (number of trials or indivoiduals under observation) n,
but exhibits local maxima and minima, depending on n, conf.level and p.hyp.
The power can be identical for different methods, depending on the particular combination of n, p.hyp, conf.level.
Note, that coverage probability and power are not necessarily symmetric for upper and lower bound of binomial CI, 
especially for Wald, Wilson Score and Agresti-Coull CI.
 
}
\value{ A list containing
  
  \item{power}{the power which is achieved for the specified parameters and method}
}

\author{
Frank Schaarschmidt
}

\seealso{ \code{\link{binDesign}} for iteration of a sample size n for which a specified power is reached  }
\examples{

# What is the probability to reject the null hypothesis
# H0: p >= 0.02 in order to show that the alternative
# hypothesis HA: p < 0.02 is very likely in the first
# example of if 200 seeds are taken from a seed lot and
# are checked for the proportion of defectives.
# Assume a true proportion under the alternative:
# p = 0.01, i.e. a absolute difference delta = 0.01 
# to the threshold proportion p.hyp=0.02.
# The null hypothesis can be rejected if the threshold
# p.hyp=0.02 is excluded by an 95 percent upper bound of the 
# Clopper-Pearson CI. What is the probability  of this event?

binPower(n=200, delta=0.01, p.hyp=0.02,
 alternative="less", method="CP")

# Assuming a lower true proportion (if one is satisfied
# also with the situation that we can only reject H0
# in case that the seed lot has a very high purity, e.g.
# only a proportion of 0.001 defectives )

binPower(n=200, delta=0.019, p.hyp=0.02,
 alternative="less", method="CP")

# Or use a higher sample size:

binPower(n=600, delta=0.01, p.hyp=0.02,
 alternative="less", method="CP")

}

\keyword{ htest }
\concept{ sample size }
\concept{ single proportion}
\concept{ binomial }
\concept{ power }
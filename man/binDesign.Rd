\name{binDesign}
\alias{binDesign}

\title{Sample Size Iteration for One Parameter Binomial Problem}

\description{
This function increases the sample size until a maximal sample size or a prespecified power is achieved.
}
\usage{
binDesign(nmax, delta, p.hyp, conf.level = 0.95,
 power = 0.8, method = "CP", alternative = "two.sided")
}
\arguments{
  \item{nmax}{ integer, maximal number of trials (individuals under observation) allowed in the iteration }
  \item{delta}{ absolute difference between the threshold and the true proportion which shall be detectable with the specified power }
  \item{p.hyp}{ threshold proportion to test against in the hypothesis, specify as a value between 0 and 1 }
  \item{conf.level}{ Confidence level of the decision, default is 0.95 }
  \item{power}{ Level of power to be achieved to be specified as a probability between 0 and 1 }
  \item{method}{ character string specifying the CImethod to be used for evaluation, see argument method in bgtCI }

  \item{alternative}{ character string, defining the alternative hypothesis, either \kbd{'two.sided'}, \kbd{'less'} or \kbd{'greater'}
	where \kbd{'less'} calculates the probability that p.hyp is excluded by an upper confidence limit given that the true proportion is p.hyp-delta,
	\kbd{'greater'} calculates the probability that p.hyp is excluded by a lower confidence limit given that the true proportion is p.hyp+delta.
	\kbd{'two.sided'} calculates min(power(p.hyp-delta, p.hyp+delta)) for a two-sided CI, thus can result in much lower power.
 }
  }
\details{
The power of a confidence interval here is defined as the probability that a confidence interval or limit excludes the threshold parameter (\kbd{p.hyp}) of the hypothesis. 

This function increases the number of trials (number of individuals under observation) until a pre-specified power is reached.
Since the power does not increase monotone with increasing n for binomial proportions but oscillates between local maxima and minima, 
the simple iteration given here will generally result in selecting those n, for which the given CI method shows a local minimum of 
coverage if the null hypothesis is true. The power can be identical for different methods, depending on the particular combination of n, p.hyp, conf.level.

Especially for large n, the calculation time may become large (particularly for Blaker). Then only the range of sample size which is of interest
can be specified in nmax, f.e. as: nmax=c(150,300). Alternatively, the function binPower might be used instead
to calculate power and bias only for some particular combinations of n, delta, p.hyp,... . 

Note that coverage probability and power are not necessarily symmetric for upper and lower bound of binomial CI.

The results can be visualized by application of the function plot() to the returned object of this function.
}
\value{a list 
  \item{powerout}{ power value, either the power of the first n exceeding the pre-specified power,
 or the maximal power achieved in the specified range of n if the specified power has not been reached}
  \item{nout}{corresponding sample size n (the number of trials) for which the prespecified power is reached,
 or the sample size n for which the maximal power has been reached within the specified range,
 if pre-specified power has not been reached.}

and additional input and iteration values needed only for the function plot.binDesign. 
}
\references{ 

\emph{Schaarschmidt, F. (2007).} Experimental design for one-sided confidence intervals or hypothesis tests in binomial group testing. \emph{Communications in Biometry and Crop Science 2 (1), 32-40.}
 http://agrobiol.sggw.waw.pl/cbcs/
  }

\author{
Frank Schaarschmidt
}

\seealso{ \code{\link{binPower}} for calculation of power, \code{\link{plot.binDesign}} for plot of the results }
\examples{

# Find a sample size for which the power to
# reject H0: p >= 0.1 in favor of HA: p < 0.1 is
# at least 0.9 (90 percent) in case that the
# true proportion is 0.04 (i.e. an absolute delta 
# of 0.06 to the threshold proportion p.hyp=0.1)
# The exact one sided Clopper-Pearson CI shall be
# used with default confidence level = 0.95.


sasi<-binDesign( nmax=200, delta=0.06,
 p.hyp=0.1, alternative="less", method="CP", power=0.9)
sasi

#### One can also plot the result:

 plot(sasi)


# For larger sample sizes iteration can be very time consuming.
# Better to use only a smaller range of n then:

sasi<-binDesign( nmax=c(200,300), delta=0.03, p.hyp=0.1,
 alternative="less", method="CP", power=0.9)
sasi


}
\keyword{ htest }
\concept{ binomial proportion}
\concept{ sample size }
\concept{ power}
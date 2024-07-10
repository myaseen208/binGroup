\name{nDesign}
\alias{nDesign}

\title{Iterate Sample Size in One Parameter Group Testing}
\description{
Increasing number of groups (assays, bulk samples) for a fixed group size in a binomial group testing design,
until a pre-specified power is achieved. At the same time, bias of the estimator is controlled. A hypothetical
threshold proportion p.hyp and the absolute difference delta to be detected have to be specified. 

}
\usage{
nDesign(nmax, s, delta, p.hyp, conf.level = 0.95, 
power = 0.8, alternative = "two.sided", method = "CP", biasrest = 0.05)
}

\arguments{
  \item{nmax}{either a single integer giving the maximal number of individuals n allowed in the iteration,
 or a vector of two integers giving the range of n in which power shall be iterated }
  \item{s}{integer, fixed group size (number of units per group) }
  \item{delta}{absolute difference between the threshold and the true proportion which shall be detectable with specified power }
  \item{p.hyp}{threshold proportion to test against in the hypothesis, specify as a value between 0 and 1 }
  \item{conf.level}{confidence level of the decision, default is 0.95 }
  \item{power}{level of power to be achieved to be specified as a probability between 0 and 1 }
  \item{alternative}{character string, defining the alternative hypothesis, either \kbd{'two.sided'}, \kbd{'less'} or \kbd{'greater'}
	where \kbd{'less'} calculates the probability that p.hyp is excluded by an upper confidence limit for a true proportion \kbd{p.hyp} - \kbd{delta},
	\kbd{'greater'} calculates the probability that p.hyp is excluded by a lower confidence limit for a true proportion \kbd{p.hyp} + \kbd{delta}.
	\kbd{'two.sided'} calculates min(power(\kbd{p.hyp} - \kbd{delta}, \kbd{p.hyp} + \kbd{delta})) for a two.sided CI, thus can result in much lower power. }

  \item{method}{ character string specifying the CImethod to be used for evaluation, see argument method in bgtCI }
  \item{biasrest}{ value between 0 and 1 specifying the absolute bias maximally allowed }
}
\details{
The power of a confidence interval here is defined as the probability that a confidence interval or limit excludes the threshold parameter (p.hyp) of the hypothesis. 

This function increases the number of groups (i.e. number of observations or assays in binomial group testing) until a pre-specified power is reached 
or the maximal number of groups \kbd{nmax} (specified in the function call) is reached. 
Since the power does not increase monotone with increasing n for binomial proportions but oscillates between local maxima and minima, 
the simple iteration given here will generally result in selecting those n, for which the given CI method shows a local minimum of 
coverage if the null hypothesis is true.
Bias decreases monotone with increasing the number of groups (if other parameters are fixed)
The resulting Problems of chosing a number of groups which results in satisfactory power, are solved in the following manner:

In case that the pre-specified power can be reached within the given range of n, the smallest n will be returned for which at least this power is reached, as
well as the actual power for this n.

In case that the pre-specified power is not reached within the given value, that n is returned for which maximal power is achieved, and the corresponding value of power.

In case that biasrestriction is violated even for the largest n within the given range of n, simply that n will be returned for which power was largest in the given range.
Due to discreteness of binomial distribution, power can be zero for one-sided hypothesis over a range of n. 

The power can be identical for different methods, depending on the particular combination of n, s, p.hyp, conf.level.

Especially for large n, the calculation time may become large (particularly for Blaker). Then only the range of sample size which is of interest
can be specified in nmax, f.e. as: nmax=c(150,300). Alternatively, the function bgtPower might be used instead
to calculate power and bias only for some particular combinations of n, s, delta, p.hyp,... . 

}
\value{
    A list containing
  \item{nout }{the number of groups (assays or bulk samples) necessary reach the power with the specified parameters}
  \item{powerout }{the power for the specified parameters and selected number of groups n}
  \item{biasout }{the bias for the specified parameters and the selected number of groups n}
 and a number of values specified in the function call or produced in the iteration, which are only necessary to apply
the function plot() on objects of class 'nDesign' 
}
\references{
\emph{Schaarschmidt F (2007).} Experimental design for one-sided confidence intervals or hypothesis tests in binomial group testing. Communications in Biometry and Crop Science 2 (1), 32-40.
 http://agrobiol.sggw.waw.pl/cbcs/

\emph{Swallow WH (1985).} Group testing for estimating infection rates and probabilities of disease transmission. Phytopathology Vol.75, N.8, 882-889.
 }

\author{
Frank Schaarschmidt
}

\seealso{
\code{\link{plot.nDesign}} to plot the iteration of this function

\code{\link{bgtPower}}: calculation of power and bias depending on n, s, delta, p.hyp, conf.level, method
\code{\link{sDesign}}:  function for stepwise increasing group size s for a given n in order to achieve sufficient power within a biasrestriction
\code{\link{estDesign}}:  function to choose group size s according to the minimal mse of the estimator, as given in Swallow (1985) }
\examples{

## Assume one aims to show that a proportion is smaller
## 0.005 (i.e. 0.5 per cent) with a power 
## of 0.80 (i.e. 80 per cent) if the unknown proportion
## in the population is 0.003 (i.e. 0.3 per cent),
## thus, to detect a delta of 0.002.
## The Clopper Pearson CI shall be used. 
## The maximal group size because of limited
## sensitivity of assay might be s=20 and we
## can only afford to perform maximally 100 assays:

nDesign(nmax=100, s=20, delta=0.002, p.hyp=0.005,
 alternative="less", method="CP", power=0.8)

## A power of 80 per cent can not be reached but
## only 30 percent with n=100 
## One might accept to show significance only for a 
## lower true proportion = 0.001 i.e accepting to be
## able to show significance only if true proportion 
## is delta=0.004 smaller than the threshold

nDesign(nmax=100, s=20, delta=0.004, p.hyp=0.005,
 alternative="less", method="CP", power=0.8)


test<-nDesign(nmax=100, s=30, delta=0.004, p.hyp=0.005,
 alternative="less", method="CP", power=0.8)

plot(test)

}
\keyword{ htest }
\concept{ group testing}
\concept{ sample size }
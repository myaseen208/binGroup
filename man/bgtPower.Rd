\name{bgtPower}
\alias{bgtPower}
\alias{bgtPowerI}

\title{Power to Reject a Hypothesis in Binomial Group Testing for One Proportion}

\description{Closed calculation of the Power to reject a hypothesis in a binomial group testing 
experiment using confidence intervals for decision. Closed calculation of bias of the point estimator 
for a given experimental design n, s and the true, unknown proportion. }

\usage{
bgtPower(n, s, delta, p.hyp, conf.level = 0.95,
 method = "CP", alternative = "two.sided")
}
\arguments{
  \item{n}{integer, giving the number of groups i.e. assays i.e. observations, a vector of integers is also allowed }
  \item{s}{integer, giving the common size of groups i.e. the number of individual units in each group, a vector of integers is also allowed  }
  \item{delta}{absolute difference between the threshold and the true proportion which shall be detectable with specified power, a vector is also allowed  }
  \item{p.hyp}{number between 0 and 1, specifying the threshold proportion in the hypotheses }
  \item{conf.level}{confidence level required for the decision on the hypotheses }
  \item{method}{character string, specifying the Confidence interval method (see \code{\link{bgtCI}}) to be used }
  \item{alternative}{character string, defining the alternative hypothesis, either 'two.sided', 'less' or 'greater'
where 'less' calculates the probability that \kbd{p.hyp} is excluded by an upper confidence limit for a true proportion p.hyp-delta,
'greater' calculates the probability that \kbd{p.hyp} is excluded by a lower confidence limit for a true proportion p.hyp+delta,
'two.sided' calculates min(power(p.hyp-delta, p.hyp+delta)) for a two-sided CI, thus can result in much lower power.
}
}
\details{
The power of a confidence interval here is defined as the probability that a confidence interval or limit excludes the threshold parameter (p.hyp) of the null hypothesis, as described in Schaarschmidt(2007). 
I.e., the null hypothesis H0: p >= p.hyp might be rejected, if an upper confidence limit for p does not contain p.hyp.
Due to discreteness, the power does not increase monotone for increasing  number of groups n or group size s, but exhibits local maxima and minima, depending on n,s, conf.level, p.hyp.
The power can be identical for different methods, depending on the particular combination of n, s, p.hyp, conf.level.
Note that coverage probability and power are not necessarily symmetric for upper and lower bounds of binomial CI, 
especially for Wald, Wilson Score and Agresti-Coull CI.

Additional to the power, bias of the point estimator is calculated according to Swallow (1985).

If vectors are specified for n, s, and (or) delta, a matrix will be constructed and power and bias are calculated for each line in this matrix. 
}
\value{
A matrix containing the columns
  \item{ns}{a vector of total sample size n*s resulting from the latter}
  \item{n}{a vector of number of groups}
  \item{s}{a vector of group sizes}
  \item{delta}{a vector of delta}
  \item{power}{ the power to reject the given null hypothesis, with the specified method and parameters of the first 4 columns}
  \item{bias}{the bias of the estimator for the specified n, s, and the true proportion}
}
\references{
\emph{Schaarschmidt F (2007).} Experimental design for one-sided confidence intervals or hypothesis tests in binomial group testing. Communications in Biometry and Crop Science 2 (1), 32-40.
 http://agrobiol.sggw.waw.pl/cbcs/

\emph{Swallow WH (1985).} Group testing for estimating infection rates and probabilities of disease transmission. Phytopathology Vol.75, N.8, 882-889.
}


\seealso{ 
\code{\link{nDesign}}:  stepwise increasing n (for a fixed group size s) until a certain power is reached within a restriction of bias of the estimator
\code{\link{sDesign}}:  stepwise increasing s (for a fixed number of groups) until a certain power is reached within a restriction of bias of the estimator 
\code{\link{estDesign}}: selection of an appropriate design to achieve minimal mean square error of the estimator
 }
\examples{

# Calculate the power for the design
# in the example given in Tebbs and Bilder(2004):
# n=24 groups each containing 7 insects
# if the true proportion of virus vectors
# in the population would be 0.04 (4 percent),
# the power to reject H0: p>=0.1 using an
# upper Clopper-Pearson ("CP") confidence interval
# can be calculated using the following call:

bgtPower(n=24, s=7, delta=0.06, p.hyp=0.1,
 conf.level=0.95, alternative="less", method="CP")


# c(), seq() or rep() might be used to explore development
# of power and bias for varying n, s, delta. How much can
# we decrease the number of groups (costly assays to be performed)
# by pooling the same number of 320 individuals to groups of
# increasing size without largely decreasing power?

bgtPower(n=c(320,160,80,64,40,32,20,10,5),
 s=c(1,2,4,5,8,10,16,32,64),
 delta=0.01, p.hyp=0.02)


# How does power develop for increasing differences
# delta between the true proportion and the threshold proportion?

bgtPower(n=50, s=10, delta=seq(from=0, to=0.01, by=0.001),
 p.hyp=0.01, method="CP")


# use a more liberal method:

bgtPower(n=50, s=10, delta=seq(from=0, to=0.01, by=0.001),
 p.hyp=0.01, method="SOC")

}

\keyword{ htest }
\concept{ sample size }
\concept{ group testing }
\concept{ binomial }
\concept{ power }

\name{bgtWidth}
\alias{bgtWidth}
\alias{bgtWidthI}

\title{ Expected Width of Confidence Intervals in Binomial Group Testing }

\description{ Calculation of expected value of the width of confidence intervals for one proportion in binomial group testing,
in dependence of the number of groups, group size, confidence level and an assumed true proportion.
Available for the confidence interval methods in bgtCI(binGroup).  
}
\usage{
bgtWidth(n, s, p, conf.level = 0.95, alternative = "two.sided",
 method = "CP")
}
\arguments{
  \item{n}{ integer, giving the number of groups i.e. assays i.e. observations, vector of integers is also allowed}
  \item{s}{ integer, giving the common size of groups i.e. the number of individual units in each group, vector of integers is also allowed }
  \item{p}{ assumed true proportion of individuals showing the trait to be estimated, vector is also allowed }
  \item{conf.level}{ required confidence level of the interval }
  \item{alternative}{ character string, defining the alternative hypothesis, either 'two.sided', 'less' or 'greater'
where 'less' calculates the expected width between the assumed true proportion p and the upper conf.level*100 percent-bound of a one-sided CI,
'greater' calculates the expected width between the true assumed proportion p and the lower conf.level*100 percent-bound of a one-sided CI,
'two.sided' calculates the expected width between the lower and the upper bound of a two-sided conf.level*100 percent-CI.}
 
 \item{method}{ character string as in the method argumnet in \code{\link{bgtCI}} }

}
\details{
For calculation of expected interval width in the standard binomial estimation see, e.g., Brown et al. (2001).
The calculation in case of binomial group testing is simply done by replacing the binomial probabilities by the
probabilities P(Y=y) for group testing (see Tebbs and Bilder, 2004)  


}
\value{ A matrix containing the columns
  \item{ns}{the resulting total number of units n*s}
  \item{n}{number of groups} 
  \item{s}{group size}
  \item{p}{the assumed true proportion }
and the calculated 
  \item{expCIWidth}{expected value of CI width as defined under argument alternative}
}
\references{
\emph{Brown LD, Cai TT, DasGupta A (2001)} Interval estimation for a binomial proportion. \emph{Statistical Science 16 (2), 101-133.}

\emph{Schaarschmidt F (2007)} Experimental design for one-sided confidence intervals or hypothesis tests in binomial group testing. \emph{Communications in Biometry and Crop Science 2 (1), 32-40.}
 http://agrobiol.sggw.waw.pl/cbcs/

\emph{Tebbs JM & Bilder CR (2004)} Confidence interval procedures for the probability of disease transmission in multiple-vector-transfer designs. \emph{Journal of Agricultural, Biological and Environmental Statistics 9 (1), 75-90.}


 }

\author{
Frank Schaarschmidt
}

\examples{

# There is a minimal expected CI length, if 
# group size s is increased (fixed other parameters)
# the corresponding group size might be chosen:

bgtWidth(n=20, s=seq(from=1, to=200, by=10),
 p=0.01, alternative="less", method="CP" )

# and this depends largely on the assumed proportion p: 

bgtWidth(n=20, s=seq(from=1, to=200, by=10),
 p=0.05, alternative="less", method="CP" )

bgtWidth(n=20, s=seq(from=1, to=200, by=10),
 p=0.005, alternative="less", method="CP" )
}

\keyword{ htest }
\concept{ sample size }
\concept{ group testing }
\concept{ binomial }
\concept{ interval width }

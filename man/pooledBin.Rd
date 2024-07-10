\name{pooledBin}
\alias{pooledBin}

\title{
Confidence intervals for a single proportion
}
\description{
Calculates confidence intervals for a single proportion based on pooled testing experiments containing various different pool sizes.
}
\usage{
pooledBin(x, m, n = rep(1, length(x)),
 pt.method = c("firth", "gart", "bc-mle", "mle", "mir"),
 ci.method = c("skew-score", "bc-skew-score", "score",
 "lrt", "wald", "mir"),
 scale = 1, alpha = 0.05, tol = .Machine$double.eps^0.5)
}
\arguments{
  \item{x}{ a vector, specifying the observed number of positive pools, among the number of pools tested (\code{n})}
  \item{m}{ a vector of pool sizes, must have the same length as, \code{x}}
  \item{n}{ a vector of the corresponding number of pools of sizes \code{m}}
  \item{pt.method}{ a character string, specifying the point estimate to compute, with the following options: \code{"bc-mle"} bias-corrected MLE, the default; \code{"mle"} MLE, and \code{"mir"} MIR. }
  \item{ci.method}{ a character string, specifying the confidence interval to compute, with options: \code{"skew-score"} skewness-corrected, the default, \code{"score"} the score, \code{"bc-skew-score"} bias- and skewness-corrected 
 \code{"lrt"} likelihood ratio test, \code{"wald"} wald, and \code{"mir"} MIR.}
  \item{scale}{ a single numeric, coefficient to scale the point estimates and intervals bounds in the print and summary method (\code{\link{print.poolbin}}, \code{\link{summary.poolbin}}) }
  \item{alpha}{ a single numeric, specifying the type-I-error level}
  \item{tol}{accuracy required for iterations in internal functions}
}


\details{
Point estimation: bias preventative ("firth") and bias corrected ("gart") estimators are recommended, with details described in Hepworth G, Biggerstaff BJ (2017). Use of MLE (\code{"mle"}), and MIR (\code{"mir"}) estimators is not recommended.

Confidence intervals: Note, that the methods "mir" and "wald" can not be recommended, because they return too narrow intervals in relevant situations, "mir" because it ignores the pooling, and "wald" because it relies on simple large sample methods.
 For computational details and simulation results of the remaining methods, see Biggerstaff (2008).}


\value{
A list with elements

\item{p}{the estimated proportion}
\item{lcl}{the lower confidence limit}
\item{ucl}{the upper confidence limit}
\item{pt.method}{the method used for point estimation }
\item{ci.method}{the method used for interval estimation }
\item{alpha}{the type-I-error level }
\item{x}{the numbers of postive pools }
\item{m}{the size of the pools }
\item{n}{the numbers of pools with corresponding pool sizes m}
\item{scale}{Scaling coefficient for the output}

}
\references{

\emph{Walter SD, Hildreth SW, Beaty BJ:} Estimation of infection rates in population of organisms using pools of variable size. \emph{Am J Epidemiol 1980, 112(1):124-128}

\emph{Hepworth G: Estimation of proportions by group testing. PhD Dissertation. 
   Melbourne, Australia: The University of Melbourne; 1999.}

\emph{Biggerstaff BJ (2008)}: Confidence interval for the difference of proportions estmimated from pooled samples. \emph{Journal of Agricultural Biological and Environmental Statistics 2008, 13(4):478-496.}

\emph{Hepworth G, Biggerstaff BJ:} Bias correction in estimating Proportions by pooled testing. \emph{JABES 20178, to appear.}
}

\author{Brad Biggerstaff}

\seealso{
 \code{\link{bgtvs}} to compute exact confidence intervals for one proportion in designs with different pool sizes. Note that \code{bgtvs} can only deal with a limited number of different pool sizes.

 \code{\link{bgtCI}} to compute exact or asymptotic confidence intervals for one proportion in designs with a common pool size in all pools.
}


\examples{

# Consider an imaginary example, where pools of size
# 1, 5, 10 and 50 are tested, 5 pools of each size
# among the 5 pools with size 1 and 5, no pool is positive,
# while among the 5 pools of size 10 and 50, 1 and 2 positive
# pools are identified, respectively. 

x1 <- c(0,0,1,2)
m1 <- c(1,5,10,50)
n1 <- c(5,5,5,5)

pooledBin(x=x1, m=m1, n=n1)
pooledBin(x=x1, m=m1, n=n1, scale=1000)

pooledBin(x=x1, m=m1, n=n1)

summary(pooledBin(x=x1, m=m1, n=n1), scale=1000)

# For another population, tested with the same design, one might find:
# 1 positive result among the pools pooling 5 elements,
# no positive result among the pools pooling 10 elements,
# 4 positive results among the pools pooling 50 elements,

x2<-c(0,1,0,4)
m2 <- c(1,5,10,50)
n2 <- c(5,5,5,5)

pooledBin(x=x2, m=m2, n=n2)

# Some other methods for the confidence bounds:

pooledBin(x=x2, m=m2, n=n2, ci.method="lrt")


###
# Reproducing some of the estimates from Table 1 in
# Hepworth & Biggerstaff (2017):

pooledBin(x=c(1,2), m=c(20,5), n=c(8,8), pt.method="firth", ci.method="lrt")

pooledBin(x=c(7,8), m=c(20,5), n=c(8,8), pt.method="firth", ci.method="lrt")

}

\keyword{ htest }
\concept{ pooled testing }
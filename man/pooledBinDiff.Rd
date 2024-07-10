\name{pooledBinDiff}
\alias{pooledBinDiff}

\title{
Confidence intervals for the difference of proportions
}
\description{
The function calculates confidence intervals for the difference of two proportions based on pooled testing experiments containing various different pool sizes.
}
\usage{
pooledBinDiff(x1, m1, x2, m2, n1 = rep(1, length(x1)),
 n2 = rep(1, length(x2)), 
 pt.method = c("firth", "gart", "bc-mle", "mle", "mir"), 
 ci.method = c("skew-score", "bc-skew-score", "score", "lrt", "wald", "mir"),
 scale = 1, alpha = 0.05, tol = .Machine$double.eps^0.5)
}

\arguments{
  \item{x1}{ a vector, specifying the observed number of positive pools, among the number of pools tested (\code{n1}) in population 1}
  \item{m1}{ a vector of corresponding pool sizes in population 1, must have the same length as \code{x1}}
  \item{x2}{ a vector, specifying the observed number of positive pools, among the number of pools tested (\code{n2}) in population 2}
  \item{m2}{ a vector of corresponding pool sizes in population 2, must have the same length as \code{x2}}
  \item{n1}{ a vector of the corresponding number of pools with sizes \code{m1}}
  \item{n2}{ a vector of the corresponding number of pools with sizes \code{m2}}
  \item{pt.method}{a character string, specifying the point estimator to compute, with the following options: \code{"firth"}: recommended bias-corrected estimator (Hepworth & Biggerstaff, 2017), the default; \code{"gart"} bias-corrected MLE (Biggerstaff, 2008); \code{"bc-mle"}: same as \code{"gart"}, kept for backward compatibility ; \code{"mle"} MLE, and \code{"mir"} MIR. }
  \item{ci.method}{a character string, specifying the confidence interval to compute, with options: \code{"skew-score"} skewness-corrected, the default, \code{"score"} the score, \code{"bc-skew-score"} bias- and skewness-corrected 
 \code{"lrt"} likelihood ratio test. See Biggerstaff (2008) for details. Further, the options \code{"wald"} wald and \code{"mir"} MIR are available but are not recommended.}
  \item{scale}{ a single numeric, coefficient to scale the point estimates and intervals bounds in the print and summary method (\code{\link{print.poolbindiff}}, \code{\link{summary.poolbindiff}}) }
  \item{alpha}{ a single numeric, specifying the type-I-error level}
  \item{tol}{accuracy required for iterations in internal functions}
}

\details{
Point estimation: bias preventative ("firth") and bias corrected ("gart") estimators are recommended, with details described in Hepworth G, Biggerstaff BJ (2017). Use of MLE (\code{"mle"}), and MIR (\code{"mir"}) estimators is not recommended.

Confidence intervals: Note, that the methods "mir" and "wald" can not be recommended, because they return too narrow intervals in relevant situations, "mir" because it ignores the pooling, and "wald" because it relies on simple large sample methods.
 For computational details and simulation results of the remaining methods, see Biggerstaff (2008).}

\value{
A list with elements

\item{d}{the estimated difference of proportions}
\item{lcl}{the lower confidence limit}
\item{ucl}{the upper confidence limit}
\item{pt.method}{the method used for point estimation }
\item{ci.method}{the method used for interval estimation }
\item{alpha}{the type-I-error level }
\item{scale}{Scaling coefficient for the output}
\item{x1}{the numbers of postive pools in population 1}
\item{m1}{the size of the pools in population 1}
\item{n1}{the numbers of pools with corresponding pool sizes m1 in population 1}
\item{x2}{the numbers of postive pools in population 2}
\item{m2}{the size of the pools in population 2}
\item{n2}{the numbers of pools with corresponding pool sizes m1 in population 2}

}
\references{

\emph{Walter SD, Hildreth SW, Beaty BJ:} Estimation of infection rates in population of organisms using pools of variable size. \emph{Am J Epidemiol 1980, 112(1):124-128}

\emph{Hepworth G: Estimation of proportions by group testing. PhD Dissertation. 
   Melbourne, Australia: The University of Melbourne; 1999.}

\emph{Biggerstaff BJ (2008)}: Confidence interval for the difference of proportions estmimated from pooled samples. \emph{Journal of Agricultural Biological and Environmental Statistics 2008, 13(4):478-496.}

\emph{Hepworth G, Biggerstaff BJ:} Bias correction in estimating Proportions by pooled testing. \emph{JABES 2017, to appear.}

}
\author{Brad Biggerstaff}

\examples{

# Consider an imaginary example, where pools of size
# 1, 5, 10 and 50 are tested, 5 pools of each size. The same 
# design is applied to estimate the prevalence in
# two populations: 
# In population 1, among the 5 pools with size 1 and 5,
# no positive pool is observed,
# while among the 5 pools of size 10 and 50, 1 and 2 positive
# pools are identified, respectively.
 
# In population 2, 1 positive result is observed
# among the 5 pools each pooling 5 elements,
# no positive result among the pools pooling 10 elements,
# 4 positive results among the pools pooling 50 elements,

x1 <- c(0,0,1,2)
m <- c(1,5,10,50)
n <- c(5,5,5,5)

x2<-c(0,1,0,4)

pooledBinDiff(x1=x1, m1=m, x2=x2, m2=m, n1=n, n2=n)

summary(pooledBinDiff(x1=x1, m1=m, x2=x2, m2= m, n1=n, n2=n), scale=1000)

# Compare recommended methods:

pooledBinDiff(x1=x1, m1=m, x2=x2, m2= m, n1=n, n2=n, pt.method="mle",
 ci.method="lrt")

pooledBinDiff(x1=x1, m1=m, x2=x2, m2= m, n1=n, n2=n, pt.method="mle",
 ci.method="score") 

pooledBinDiff(x1=x1, m1=m, x2=x2, m2= m, n1=n, n2=n, pt.method="mle",
 ci.method="skew-score")
  
}

\keyword{ htest }
\keyword{ pooled testing }

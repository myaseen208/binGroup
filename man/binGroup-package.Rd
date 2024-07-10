\name{binGroup-package}
\Rdversion{2.2-1}
\alias{binGroup-package}
\alias{binGroup}
\docType{package}
\title{
Statistical Methods for Group Testing.
}
\description{
This package provides methods for estimation and hypothesis testing of proportions in group testing designs. It involves methods for estimating a proportion in a single population (assuming sensitivity and specificity 1 in designs with equal group sizes), as well as hypothesis tests and functions for experimental design for this situation. For estimating one proportion or the difference of proportions, a number of confidence interval methods are included, which can deal with various different pool sizes. Further, regression methods are implemented for simple pooling and matrix pooling designs. \cr \cr
This package also provides methods for identification of positive items in a number of group testing algorithms. Optimal testing configurations can be found for hierarchical and array-based algorithms. Operating characteristics can be calculated for testing configurations across a wide variety of situations.
}
\details{
\tabular{ll}{
Package: \tab binGroup\cr
Type: \tab Package\cr
Version: \tab 2.2-1\cr
Date: \tab 2018-07-07\cr
License: \tab GPL (>=3)\cr
LazyLoad: \tab no\cr
}

1) One-sample case

Methods for calculating confidence intervals for a single population proportion from designs with equal group sizes (as described by Tebbs and Bilder, 2004 and Schaarschmidt, 2007) are implemented in the function \code{\link{bgtCI}}.

For the problem of choosing an adequate experimental design in the one-sample case with only one group size, the functions \code{\link{estDesign}}, \code{\link{sDesign}}, \code{\link{nDesign}} implement different iterative approaches, as examplified by Swallow (1985) and Schaarschmidt (2007).

If a confidence interval for a single proportion shall be calculated based on a design involving groups of different group sizes, a number of methods described by Hepworth (1999) is available in the function \code{\link{pooledBin}}. The exact method described by Hepworth (1996) is implemented in the function \code{\link{bgtvs}}.

2) Two-sample case

The function \code{\link{pooledBinDiff}} provides a number of confidence interval methods for estimating the difference of proportions from two independent samples, allowing for groups of different group size (Biggerstaff, 2008).

3) Regression models

Two approaches (by Vansteelandt et al., 2000 and Xie, 2001) to estimate parameters of group testing regression models can be applied by calling \code{\link{gtreg}}. Once fitted, corresponding methods to extract residuals, calculate predictions and summarize the parameter estimates (including hypotheses tests) are available in the S3 methods \code{\link{residuals.gt}}, \code{\link{predict.gt}} and \code{\link{summary.gt}}.

Group testing regression models in settings with matrix pooling (Xie, 2001) can be fit using \code{\link{gtreg.mp}}. 

4) Identification using hierarchical and array-based group testing algorithms

The function \code{\link{OTC}} implements a number of group testing algorithms, described in Hitt et al. (2018), which calculate the operating characteristics and find the optimal testing configuration over a range of possible initial pool sizes and/or testing configurations (sets of subsequent pool sizes).



}
\author{
Boan Zhang, Christopher Bilder, Brad Biggerstaff, Brianna Hitt, Frank Schaarschmidt

Maintainer: Frank Schaarschmidt <schaarschmidt@biostat.uni-hannover.de>
}

\references{
\emph{Biggerstaff, B.J. (2008)}: Confidence interval for the difference of proportions estmimated from pooled samples. \emph{Journal of Agricultural Biological and Environmental Statistics, 13(4), 478-496.}

\emph{Hepworth, G. (1996)} Exact confidence intervals for proportions estimated by group testing. \emph{Biometrics 52, 1134-1146.}

\emph{Hepworth, G. (1999)}: \emph{Estimation of proportions by group testing.} PhD Dissertation. Melbourne, Australia: The University of Melbourne.

\emph{Hitt, B., Bilder, C., Tebbs, J. & McMahan, C. (2018)} The Optimal Group Size Controversy for Infectious Disease Testing: Much Ado About Nothing?!. Manuscript submitted for publication. http://www.chrisbilder.com/grouptesting

\emph{Schaarschmidt, F. (2007)} Experimental design for one-sided confidence intervals or hypothesis tests in binomial group testing.  \emph{Communications in Biometry and Crop Science 2 (1), 32-40.} http://agrobiol.sggw.waw.pl/cbcs/

\emph{Swallow, W.H. (1985)} Group testing for estimating infection rates and probabilities of disease transmission.  \emph{Phytopathology 75 (8), 882-889.}

\emph{Tebbs, J.M. & Bilder, C.R. (2004)} Confidence interval procedures for the probability of disease transmission in multiple-vector-transfer designs.  \emph{Journal of Agricultural, Biological and Environmental Statistics 9 (1), 75-90.}

\emph{Vansteelandt, S., Goetghebeur, E., and Verstraeten, T. (2000)} Regression models for disease prevalence with diagnostic tests on pools of serum samples, \emph{Biometrics, 56, 1126-1133.}

\emph{Xie, M. (2001)} Regression analysis of group testing samples, \emph{Statistics in Medicine, 20, 1957-1969.}

}

\keyword{package}
\keyword{regression}
\keyword{models}
\keyword{htest}
\keyword{Categorical Data}
\keyword{Optimization}

\examples{

# 1) One-sample problem

# 1.1) Confidence intervals for designs with equal group size (pool size),
# where 
#  n denotes the number of groups (pools),
#  s denotes the common group size (number of individuals pooled per group),
#  y denotes the number of groups tested positive.

# The following call reproduces the example given 
# by Tebbs and Bilder (2004) for the two-sided 95-percent 
# exact (Clopper-Pearson) interval:

bgtCI(n=24, y=3, s=7, conf.level=0.95,
 alternative="two.sided", method="CP")

# 1.2) Confidence intervals for designs with unequal group size (pool size):
# Keeping notation as above but allowing for (a limited number of) different
# group size s, the examples given in Hepworth (1996), Table 5 can be 
# reproduced by calling:

 bgtvs(n=c(2,3), s=c(5,2), y=c(0,0))
 bgtvs(n=c(2,3), s=c(5,2), y=c(0,1))

# The function pooledBin provides different methods for the same problem,
# where x is the number of positive groups, m is the size of the groups and
# n is the number of groups with the correesponding sizes:

pooledBin(x=c(0,1), m=c(5,2), n=c(2,3), ci.method="score")
pooledBin(x=c(0,1), m=c(5,2), n=c(2,3), ci.method="lrt")
pooledBin(x=c(0,1), m=c(5,2), n=c(2,3), ci.method="bc-skew-score")

# 1.3) For experimental design based on the bias of the point estimate,
# as proposed by Swallow (1985): The values in Table 1 (Swallow, 1985),
# p.885 can be reproduced by calling:

estDesign(n=10, smax=100, p.tr=0.001)
estDesign(n=10, smax=100, p.tr=0.01)

# 2) Two-sample comparison

# Assume a design, where pools 5, 1, 1, 30, and 20 pools of size 10, 4, 1, 25, 50, 
# respectively, are used to estimate the prevalence in two populations.
# In population 1, one out of 5 pools with 10 units is positive,
# while in population 2, two out of five pools with 10 units is positive as well as 
# the one pool with only 1 unit.
# The difference of proportions is to be estimated.

x1 <- c(1,0,0,0,0)
m1 <- c(10,4,1,25,50)
n1 <- c(5,1,1,30,20)

x2 <- c(2,0,1,0,0)
m2 <- c(10,4,1,25,50)
n2 <- c(5,1,1,30,20)

pooledBinDiff(x1=x1, m1=m1,x2=x2, m2=m2, n1=n1, n2=n2, ci.method="lrt")


# 3) Regression models

# 3.1) Fitting a regression model
# A HIV surveillance data (used by Vansteelandt et al. 2000)
# can be analysed for the dependence of HIV prevalence
# on covariates AGE and EDUC., with sensitivity and specificity
# assumed to be 0.9 each.

data(hivsurv)
fit1 <- gtreg(formula = groupres ~ AGE + EDUC., data = hivsurv,
  groupn = gnum, sens = 0.9, spec = 0.9, method = "Xie")
summary(fit1)


# 3.2) Fitting a regression model for matrix pooling data
# The function sim.mp is used to simulate a matrix pooling data set:

set.seed(9128)
sa1a<-sim.mp(par=c(-7,0.1), n.row=c(5,4), n.col=c(6,5),
 sens=0.95, spec=0.95)

str(sa1a)
sa1<-sa1a$dframe


\dontrun{
fit2 <- gtreg.mp(formula = cbind(col.resp, row.resp) ~ x, data = sa1, 
                 coln = coln, rown = rown, arrayn = arrayn, 
                 sens = 0.95, spec = 0.95, n.gibbs = 2000, trace = TRUE)
fit2
summary(fit2)
}

# 4) Identification using hierarchical and array-based group testing algorithms

# 4.1) Finding the optimal testing configuration over a range of initial 
# group sizes, using non-informative three-stage hierarchical testing, where 
#   p denotes the overall prevalence of disease,
#   Se denotes the sensitivity of the diagnostic test, 
#   Sp denotes the specificity of the diagnostic test,
#   group.sz denotes the range of initial pool sizes for consideration, and
#   obj.fn specifies the objective functions for which to find results.

# The following call reproduces results given by Hitt et al. (2018) for 
# informative three-stage hierarchical testing with an overall disease 
# prevalence of E(p_i) = 0.01 and sensitivity and specificity equal to 0.95.

# This example takes approximately 2.5 minutes to run.

\dontrun{
set.seed(1002)
results1 <- OTC(algorithm="ID3", p=0.01, Se=0.95, Sp=0.95, group.sz=3:40, 
obj.fn=c("ET", "MAR"), alpha=2)
results1$opt.ET$OTC
results1$opt.ET$ET/results1$opt.ET$OTC$Stage1
results1$opt.MAR$OTC
results1$opt.MAR$ET/results1$opt.MAR$OTC$Stage1
}

# 4.2) Finding the optimal testing configuration using non-informative
# array testing without master pooling

# The following call reproduces results given by Hitt et al. (2018) for
# non-informative array testing without master pooling with an overall
# disease prevalence of p=0.01 and sensitivity and specificity equal 
# to 0.90.

# This example takes approximately 7 minutes to run.

\dontrun{
results2 <- OTC(algorithm="A2", p=0.01, Se=0.90, Sp=0.90, group.sz=3:40, 
obj.fn=c("ET", "MAR"))
results2$opt.ET$OTC
results2$opt.ET$ET/results2$opt.ET$OTC$Array.sz
results2$opt.MAR$OTC
results2$opt.MAR$ET/results2$opt.MAR$OTC$Array.sz
}

}

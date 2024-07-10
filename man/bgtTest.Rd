\name{bgtTest}
\alias{bgtTest}


\title{ Hypothesis Test for One Proportion in Binomial Group Testing }

\description{ Calculates p values for hypotheses tests of binomial proportions 
estimated from binomial group testing experiments against a threshold proportion
in the hypotheses. Exact test, Score test and Wald test are available methods. 
Assumes equal group sizes, 100 percent sensitivity and specificity of the assays to test the groups,
and individuals units randomly assigned to the groups with identical true probability of success.
 }

\usage{
bgtTest(n, y, s, p.hyp, alternative = "two.sided",
 method = "Exact")
}

\arguments{
  \item{n}{integer, number of  groups (i.e. assays i.e. observations)}
  \item{y}{integer, number of positive groups }
  \item{s}{integer, common size of groups i.e. the number of individual units in each group}
  \item{p.hyp}{number between 0 and 1, specifying the hypothetical threshold proportion to test against}
  \item{alternative}{ character string defining the alternative hypothesis, either 'two.sided', 'less' or 'greater' }
  \item{method}{ character string defining the test method to be used:
can be one of 'Exact' for an exact test corresponding to the Clopper-Pearson confidence interval
'Score' for a Score test, corresponding to the Wilson confidence interval
'Wald' for a Wald test corresponding to the Wald confidence interval, not recommended}
}


\value{A list containing:
  \item{p.value}{the p value of the test}
  \item{estimate}{the estimated proportion}
  \item{p.hyp}{as input}
  \item{alternative}{as input}
  \item{method}{as input}
}
\references{ 

\emph{Swallow WH (1985)} Group testing for estimating infection rates and probabilities of disease transmission. \emph{Phytopathology 75 (8), 882-889.}

\emph{Blyth, C and Still, H. (1983)} Binomial confidence intervals. \emph{Journal of the American Statistical Association 78, 108-116.}

\emph{Santner TJ and Duffy DE (1989)} The statistical analysis of discrete data. Springer New York.

\emph{Remund KM, Dixon DA, Wright DL, Holden LR (2001)} Statistical considerations on seed purity testing on transgenic traits. \emph{Seed Science Research (11), 101-119.}

 }

\seealso{ \code{\link{bgtCI}} for confidence intervals in binomial group testing}
\examples{

# Assume the experiment: Assays are performed on
# n=10 groups, each group is a bulk sample
# of s=100 individuals, aim is to show that 
# less than 0.5 percent ('p<0.005') of the units
# of the population show a detrimental trait (positive assay).
# The assay is senstive to show a positive result if only 1
# unit in the bulk sample of 100 units is positive.
# y=1 positive assay and 9 negative assays are observed.

bgtTest(n=10,y=1,s=100,alternative="less",method="Exact",p.hyp=0.005)

# The exact test corresponds to the 
# limits of the Clopper-Pearson confidence interval
# in the example of Tebbs and Bilder(2004):

bgtTest(n=24, y=3, s=7, alternative="two.sided",
 method="Exact", p.hyp=0.0543)

bgtTest(n=24, y=3, s=7, alternative="two.sided",
 method="Exact", p.hyp=0.0038)

# Further methods:

bgtTest(n=24, y=3, s=7, alternative="two.sided",
 method="Score", p.hyp=0.0516)

bgtTest(n=24, y=3, s=7, alternative="two.sided",
 method="Wald", p.hyp=0.0401)
}

\keyword{ htest }
\concept{ binomial }
\concept{ group testing }
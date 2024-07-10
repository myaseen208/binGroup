\name{summary.gt.mp}
\Rdversion{1.1}
\alias{summary.gt.mp}

\title{
Summary Method for Group Testing Model (Matrix Pooling) Fits
}
\description{
Produce a summary list for objects of class \kbd{"gt.mp"} returned by \code{\link{gtreg.mp}}.
}
\usage{
\method{summary}{gt.mp}(object, ...)
}

\arguments{
  \item{object}{
a fitted object of class \kbd{"gt.mp"}.  
}
   \item{\dots}{currently not used.}
}
\details{
\kbd{print.summary.gt.mp} is the print function that formats the coefficients, standard errors, etc. and additionally gives 'significance stars' if \kbd{signif.stars} is \kbd{TRUE}. The \kbd{coefficients} component of the result gives the estimated coefficients and their estimated standard errors, together with their ratio. This third column is labeled z ratio using Wald tests. A fourth column gives the two-tailed p-value corresponding to the \kbd{z ratio} based on a Wald test. 
}

\value{
\kbd{summary.gt.mp} returns an object of class \kbd{"summary.gt.mp"}, a list with components 
 \item{call}{the component from \kbd{object}.}
 \item{link}{the component from \kbd{object}.}
 \item{coefficients}{the matrix of coefficients, standard errors, z-values and p-values.}
 \item{counts}{the component from \kbd{object}.}
 \item{Gibbs.sample.size}{the component from \kbd{object}.}
 \item{cov.mat}{the estimated covariance matrix of the estimated coefficients.}
}


\author{
Boan Zhang
}


\seealso{
\code{\link{gtreg.mp}} for creating an object of class \kbd{"gtreg.mp"},
\code{\link{print.summary.gt.mp}} gives some hints how to change the print-out of \kbd{summary.gt.mp}
}
\examples{

## --- Continuing the Example from  '?sim.mp' and '?gtreg.mp':
# 5*6 and 4*5 matrix

set.seed(9128)
sa1a<-sim.mp(par=c(-7,0.1), n.row=c(5,4), n.col=c(6,5),
  sens=0.95, spec=0.95)
sa1<-sa1a$dframe

\dontrun{
fit1mp <- gtreg.mp(formula = cbind(col.resp, row.resp) ~ x, data = sa1, 
  coln = coln, rown = rown, arrayn = arrayn, 
  sens = 0.95, spec = 0.95, linkf = "logit", n.gibbs = 1000, tol = 0.005)

summary(fit1mp)
}

}

\keyword{print}

\name{summary.gt}
\Rdversion{1.1}
\alias{summary.gt}

\title{
Summary Method for Group Testing Model (Simple Pooling) Fits
}
\description{
Produce a summary list for objects of class \kbd{"gt"} returned by \kbd{gtreg}.
}
\usage{
\method{summary}{gt}(object, ...)
}

\arguments{
  \item{object}{
a fitted object of class \kbd{"gt"}.  
}
   \item{\dots}{currently not used.}
}
\details{
\kbd{print.summary.gt} is the print function that formats the coefficients, standard errors, etc. and additionally gives 'significance stars' if \kbd{signif.stars} is \kbd{TRUE}. The \kbd{coefficients} component of the result gives the estimated coefficients and their estimated standard errors, together with their ratio. This third column is labeled \kbd{z ratio} using Wald tests. A fourth column gives the two-tailed p-value corresponding to the z ratio based on a Wald test. (It is possible that there are no residual degrees of freedom from which to estimate it. In that case the estimate is \kbd{NaN}.) 
}
\value{
\kbd{summary.gt} returns an object of class \kbd{"summary.gt"}, a list with components 
 \item{call}{the component from \kbd{object}.}
 \item{link}{the component from \kbd{object}.}
 \item{deviance}{the component from \kbd{object}.}
 \item{aic}{the component from \kbd{object}.}
 \item{df.residual}{the component from \kbd{object}.}
 \item{null.deviance}{the component from \kbd{object}.}
 \item{df.null}{the component from \kbd{object}.}
 \item{deviance.resid}{the deviance residuals: see \code{\link{residuals.gt}}}
 \item{coefficients}{the matrix of coefficients, standard errors, z-values and p-values. Aliased coefficients are omitted.}
 \item{counts}{the component from \kbd{object}.}
 \item{method}{the component from \kbd{object}.}
 \item{cov.mat}{the estimated covariance matrix of the estimated coefficients.}
}

\author{
Boan Zhang
}

\seealso{
 \code{\link{gtreg}} for creating an object of class \kbd{"gt"}, 
  and \code{\link{print.summary.gt}} for some options of changing the print out following \kbd{summary.gt}.
}
\examples{
## --- Continuing the Example from '?gtreg':

data(hivsurv)

fit1 <- gtreg(formula = groupres ~ AGE + EDUC., data = hivsurv,
  groupn = gnum, sens = 0.9, spec = 0.9, method = "Xie")

summary(fit1)

}

\keyword{print}
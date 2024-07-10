\name{residuals.gt}
\Rdversion{1.1}
\alias{residuals.gt}

\title{
Extract Model Residuals From a Fitted Group Testing Model
}
\description{
\kbd{residuals.gt} is a function which extracts model residuals from objects of class \kbd{"gt"} returned by \code{\link{gtreg}}. 
}
\usage{
\method{residuals}{gt}(object, type = c("deviance", "pearson", "response"),...)
}
\arguments{
  \item{object}{
an object of class \kbd{"gt"} from which the model residuals are to be extracted.
}
  \item{type}{
the type of residuals which should be returned. The alternatives are: \kbd{"deviance"} (default), \kbd{"pearson"}, and \kbd{"response"}.
}
\item{\dots}{currently not used}
}

\value{
Residuals of group responses extracted from the object \kbd{object}.
}

\author{
Boan Zhang
}

\examples{

data(hivsurv)

fit1 <- gtreg(formula = groupres ~ AGE * EDUC., data = hivsurv, groupn = gnum,  
              linkf = "probit")
residuals.gt(object = fit1, type = "pearson")
residuals.gt(object = fit1, type = "deviance")

}

\keyword{regression}
\keyword{models}
\concept{residuals}

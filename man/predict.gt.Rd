\name{predict.gt}
\Rdversion{1.1}
\alias{predict.gt}

\title{
Predict Method for Group Testing Model Fits
}
\description{
Obtains predictions for individual observations and optionally estimates standard errors of those predictions from objects of class \kbd{"gt"} or \kbd{"gt.mp"} returned by \kbd{gtreg} and \kbd{gtreg.mp}, respectively.
}
\usage{
\method{predict}{gt}(object, newdata, type = c("link", "response"),
 se.fit = FALSE, conf.level = NULL, na.action = na.pass, ...)
}

\arguments{
  \item{object}{
a fitted object of class \kbd{"gt"} or \kbd{"gt.mp"}.  
}
  \item{newdata}{
optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.
}
  \item{type}{
the type of prediction required. The option \kbd{"link"} is on the scale of the linear predictors; the alternative \kbd{"response"} is on the scale of the response variable. Thus for the binomial model the \kbd{"link"} predictions are of log-odds (probabilities on logit scale) and \kbd{type = "response"} gives the predicted probabilities. 
}
  \item{se.fit}{
logical switch indicating if standard errors are required.
}
  \item{conf.level}{
confidence level of the interval of the predicted values.
}
  \item{na.action}{
function determining what should be done with missing values in \kbd{newdata}. The default is to predict \kbd{NA}.
}
  \item{\dots}{ currently not used}
}
\details{
If \kbd{newdata} is omitted the predictions are based on the data used for the fit. When \kbd{newdata} is present and contains missing values, how the missing values will be dealt with is determined by the \kbd{na.action} argument. In this case, if \kbd{na.action = na.omit} omitted cases will not appear, whereas if \kbd{na.action = na.exclude} they will appear (in predictions and standard errors), with value \kbd{NA}. See also \kbd{napredict}. 
}
\value{
If \kbd{se = FALSE}, a vector or matrix of predictions. If \kbd{se = TRUE}, a list with components

  \item{fit}{ Predictions }
  \item{se.fit}{ Estimated standard errors }
  \item{lower}{ Lower bound of the confidence interval if calculated }
  \item{upper}{ Upper bound of the confidence interval if calculated }

}

\author{
Boan Zhang
}

\examples{

data(hivsurv)

fit1 <- gtreg(formula = groupres ~ AGE + EDUC., data = hivsurv,  
  groupn = gnum, sens = 0.9, spec = 0.9, linkf = "logit", method = "V")
pred.data <- data.frame(AGE = c(15, 25, 30), EDUC. = c(1, 3, 2))
predict(object = fit1, newdata = pred.data, type = "link", se.fit = TRUE)
predict(object = fit1, newdata = pred.data, type = "response", 
  se.fit = TRUE, conf.level = 0.9)
predict(object = fit1, type = "response", se.fit = TRUE, conf.level = 0.9)


}

\keyword{regression}
\keyword{models}
\concept{prediction}
\concept{group testing}
\name{gtreg}
\Rdversion{1.1}
\alias{gtreg}
\alias{gtreg.fit}
\alias{EM}
\alias{EM.ret}

\title{
Fitting Group Testing Models
}

\description{
\kbd{gtreg} is a function to fit the group testing regression model specified through a symbolic description of the linear predictor and descriptions of the group testing setting.
}

\usage{
gtreg(formula, data, groupn, retest = NULL, sens = 1,
 spec = 1, linkf = c("logit", "probit", "cloglog"),
 method = c("Vansteelandt", "Xie"), sens.ind = NULL,
 spec.ind = NULL, start = NULL, control = gt.control(...), ...)

gtreg.fit(Y, X, groupn, sens, spec, linkf, start=NULL)

EM(Y, X, groupn, sens, spec, linkf, start = NULL,
 control = gt.control())

EM.ret(Y, X, groupn, ret, sens, spec, linkf, sens.ind,
 spec.ind, start = NULL, control = gt.control())

}


\arguments{
  \item{formula}{
an object of class \kbd{"formula"} (or one that can be coerced to that class): a symbolic description of 
the model to be fitted. The details of model specification are given under 'Details'.

}
  \item{data}{
an optional data frame, list or environment (or object coercible by \kbd{as.data.frame} to a data 
frame) containing the variables in the model. If not found in data, the variables are taken from 
\kbd{environment(formula)}, typically the environment from which \kbd{gtreg} is called.

}
  \item{groupn}{
a vector, list or data frame of the group numbers that designates individuals to groups.
}
  \item{retest}{
a vector, list or data frame of individual retest results for Dorfman's retesting procedure. Default value is \kbd{NULL} for no retests. See 'Details' for how to code it.
}  

  \item{sens}{
sensitivity of the test, set to be 1 by default.
}
  \item{spec}{
specificity of the test, set to be 1 by default.
}
  \item{sens.ind}{
sensitivity of the individual retests, set to be equal to \kbd{sens} if not specified otherwise.
}
  \item{spec.ind}{
specificity of the individual retests, set to be equal to \kbd{spec} if not specified otherwise.
}

  \item{linkf}{
a character string specifying one of the three link functions for a binomial model: \kbd{"logit"} 
(default) or \kbd{"probit"} or \kbd{"cloglog"}.

}
  \item{method}{
The method to fit the model, must be one of \kbd{"Vansteelandt"} (default) or \kbd{"Xie"}. The option 
\kbd{"Vansteelandt"} finds estimates by directly maximizing the likelihood function based on the group responses while the option \kbd{"Xie"} uses the EM algorithm to maximize the likelihood function in terms of the unobserved individual responses.
}

\item{start}{
starting values for the parameters in the linear predictor.
}


\item{control}{
a list of parameters for controlling the fitting process in method \kbd{"Xie"}. See the documentation for \code{\link{gt.control}} for details.
}

\item{Y}{
For \kbd{gtreg.fit}, \kbd{EM} and \kbd{EM.ret}: the vector of the group responses.
}

\item{X}{
For \kbd{gtreg.fit}, \kbd{EM} and \kbd{EM.ret}: the design matrix of the covariates.
}

\item{ret}{
For \kbd{EM.ret}: a vector containing individual retest results.
}

\item{...}{arguments to be passed by default to \code{\link{gt.control}}: see argument \kbd{control} }


}
\details{

A typical predictor has the form \kbd{groupresp ~ covariates} where response is the (numeric) group response vector and \kbd{covariates} is a series of terms which specifies a linear predictor for individual responses. Note that it is actually the unobserved individual responses, not the observed group responses, which are modeled by the covariates here. In \kbd{groupresp}, a 0 denotes a negative response and a 1 denotes a positive response, where the probability of an individual positive response is being modeled directly.  A terms specification of the form \kbd{first + second} indicates all the terms in \kbd{first} together with all the terms in \kbd{second} with duplicates removed. The \kbd{terms} in the formula will be re-ordered so that main effects come first, followed by the interactions, all second-order, all third-order and so on; to avoid this pass a terms object as the formula.

A specification of the form \kbd{first:second} indicates the set of terms obtained by taking the interactions of all terms in \kbd{first} with all terms in \kbd{second}. The specification \kbd{first*second} indicates the cross of \kbd{first} and \kbd{second}. This is the same as \kbd{first + second + first:second}. 

Three workhorse functions \kbd{gtreg.fit}, \kbd{EM} and \kbd{EM.ret}, where the first corresponds to Vansteelandt's method and the last two corresponds to Xie's method, are called by \kbd{gtreg} to carry out the model fitting. The \kbd{gtreg.fit} function uses the \kbd{optim} function with default method \kbd{"Nelder-Mead"} to maximize the likelihood function of the observed group responses. If this optimization method produces a Hessian matrix of all zero elements, the \kbd{"SANN"} method in \kbd{optim} is employed to find the coefficients and Hessian matrix. For \kbd{"SANN"} method, the number of iterations in \kbd{optim} is set to be 10000.

The \kbd{EM} and \kbd{EM.ret} function apply Xie's EM algorithm to the likelihood function written in terms of the unobserved individual responses; the functions use \kbd{glm.fit} to update the parameter estimates within each M step. The \kbd{EM} function is used when there are no retests and \kbd{EM.ret} is used when individual retests are available. Thus, within \kbd{retest}, individual observations in observed positive groups are 0 (negative) or 1 (positive); the remaining individual observations are \kbd{NA}s meaning that no retest is performed for them. Retests cannot be used with Vansteelandt's method; a warning message will be given in this case, and the individual retests will be ignored in the model fitting. There could be slight differences in the estimates between the Vansteelandt's and Xie's methods (when retests are not available) due to different convergence criteria. 

The data used here should be in the form of simple pooling - meaning that each individual appears in exactly one pool. When only the group responses are observed, the null degrees of freedom are the number of groups minus 1 and the residual degrees of freedom are the number of groups minus the number of parameters. When individual retests are observed too, it is an open research question for what the degrees of freedom and the deviance for the null model should
be; therefore the degrees of freedom and \kbd{null.deviance} will not be displayed.

For the background on the use of \kbd{optim}, see \kbd{help(optim)}.

}
\value{
\kbd{gtreg} returns an object of class \kbd{"gt"}. See later in this section. 
The function \kbd{summary} (i.e., \code{\link{summary.gt}}) can be used to obtain or print a summary of the results. 
The group testing functions \kbd{predict} (i.e., \code{\link{predict.gt}}) and \kbd{residuals} (i.e., \code{\link{residuals.gt}}) can be used to extract various useful features of the value returned by \kbd{gtreg}. 
An object of class \kbd{"gt"} is a list containing at least the following components: 

 \item{ coefficients }{ a named vector of coefficients }
 \item{ hessian }{ estimated Hessian matrix of the negative log likelihood function, serves as an estimate of the information matrix }
 \item{ residuals }{ the response residuals, difference of the observed group responses and the fitted group responses. }
 \item{ fitted.values }{ the fitted mean values of group responses. }
 \item{ deviance }{ the deviance between the fitted model and the saturated model. }
 \item{ aic }{ Akaike's An Information Criterion, minus twice the maximized log-likelihood plus twice the number of coefficients }
 \item{ null.deviance }{ The deviance for the null model, comparable with \kbd{deviance}. The null model will include only the intercept if there is one in the model. }
 \item{ counts }{ For Vansteelandt's method: the number of iterations in \kbd{optim}; For Xie's method: the number of iterations in the EM algorithm.}
 \item{ df.residual }{ the residual degrees of freedom. }
 \item{ df.null }{ the residual degrees of freedom for the null model. }
 \item{ z }{ the vector of group responses. }
 \item{ call }{ the matched call. }
 \item{ formula }{ the formula supplied. }
 \item{ terms }{ the terms object used. }
 \item{ method }{ the method (\kbd{"Vansteelandt"} or \kbd{"Xie"}) used to fit the model. }
 \item{ link }{ the link function used in the model. }
}
\references{

\emph{Xie, M. (2001)}, Regression analysis of group testing samples, \emph{Statistics in Medicine, 20, 1957-1969.}

\emph{Vansteelandt, S., Goetghebeur, E., and Verstraeten, T. (2000)}, Regression models for disease prevalence with diagnostic tests on pools of serum samples, \emph{Biometrics, 56, 1126-1133.}

}

\author{
Boan Zhang
}


\seealso{

\code{\link{summary.gt}}, \code{\link{predict.gt}} and \code{\link{residuals.gt}} for \kbd{gt} methods.
\code{\link{gtreg.mp}} for the group testing regression model in the matrix pooling setting.


}
\examples{

data(hivsurv)

fit1 <- gtreg(formula = groupres ~ AGE + EDUC., data = hivsurv,
           groupn = gnum, sens = 0.9, spec = 0.9, method = "Xie")
fit1

## --- Continuing the Example from  '?sim.gt':

set.seed(46)
gt.data <- sim.gt(par = c(-12, 0.2), sample.size = 700, group.size = 5)
fit2 <- gtreg(formula = gres ~ x, data = gt.data, groupn = groupn)
fit2

set.seed(21)
gt.data <- sim.gt(par = c(-12, 0.2), sample.size = 700, group.size = 6,
               sens = 0.95, spec = 0.95, sens.ind = 0.98, spec.ind = 0.98)
fit1 <- gtreg(formula = gres ~ x, data = gt.data, groupn = groupn, 
        retest = retest, method = "X", sens = 0.95, spec = 0.95, sens.ind = 0.98,
        spec.ind = 0.98, trace = TRUE)
summary(fit1)


}

\keyword{models}
\keyword{regression}
\concept{group testing}
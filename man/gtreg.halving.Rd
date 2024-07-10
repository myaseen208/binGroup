\name{gtreg.halving}
\Rdversion{1.1}
\alias{gtreg.halving}
\alias{EM.halving}

\title{
Fitting Group Testing Models Under the Halving Protocol
}

\description{
\kbd{gtreg.halving} is a function to fit the group testing regression model under the halving protocol specified through a symbolic description of the linear predictor and descriptions of the group testing setting.
}

\usage{
gtreg.halving(formula, data, groupn, subg, retest, sens = 1,
 spec = 1, linkf = c("logit", "probit", "cloglog"),
 sens.ind = NULL, spec.ind = NULL, start = NULL, control = gt.control(...), ...)

EM.halving(Y, X, groupn, subg, ret, sens, spec, linkf, sens.ind,
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
\kbd{environment(formula)}, typically the environment from which \kbd{gtreg.halving} is called.
}
  \item{groupn}{
a vector, list or data frame of the group numbers that designates individuals to groups.
}
  \item{subg}{
a vector, list or data frame of the group numbers that designates individuals to subgroups.
}
  \item{retest}{
a vector, list or data frame of individual retest results within a subgroup that tests positive.
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
\item{start}{
starting values for the parameters in the linear predictor.
}


\item{control}{
a list of parameters for controlling the fitting process. See the documentation for \code{\link{gt.control}} for details.
}

\item{Y}{
For \kbd{EM.halving}: the vector of the group responses.
}

\item{X}{
For \kbd{EM.halving}: the design matrix of the covariates.
}

\item{ret}{
For \kbd{EM.halving}: a vector containing individual retest results.
}

\item{...}{arguments to be passed by default to \code{\link{gt.control}}: see argument \kbd{control} }


}
\details{

A typical predictor has the form \kbd{groupresp ~ covariates} where response is the (numeric) group response vector and \kbd{covariates} is a series of terms which specifies a linear predictor for individual responses. Note that it is actually the unobserved individual responses, not the observed group responses, which are modeled by the covariates here. In \kbd{groupresp}, a 0 denotes a negative response and a 1 denotes a positive response, where the probability of an individual positive response is being modeled directly.  A terms specification of the form \kbd{first + second} indicates all the terms in \kbd{first} together with all the terms in \kbd{second} with duplicates removed. The \kbd{terms} in the formula will be re-ordered so that main effects come first, followed by the interactions, all second-order, all third-order and so on; to avoid this pass a terms object as the formula.

A specification of the form \kbd{first:second} indicates the set of terms obtained by taking the interactions of all terms in \kbd{first} with all terms in \kbd{second}. The specification \kbd{first*second} indicates the cross of \kbd{first} and \kbd{second}. This is the same as \kbd{first + second + first:second}. 

The \kbd{EM.halving} function apply Xie's EM algorithm to the likelihood function written in terms of the unobserved individual responses; the functions use \kbd{glm.fit} to update the parameter estimates within each M step. In the halving protocol, if the initial group tests positive, it is split into two subgroups. The two subgroups are subsequently tested and if either subgroup tests positive, the third and final step is to test all individuals within the subgroup. Thus, within \kbd{subg}, subgroup responses in observed positive groups are 0 (negative) or 1 (positive); the remaining subgroup responses are \kbd{NA}s meaning that no tests are performed for them. The individual retests are similarly coded. 

}
\value{
\kbd{gtreg.halving} returns an object of class \kbd{"gt"}. See later in this section. 
The function \kbd{summary} (i.e., \code{\link{summary.gt}}) can be used to obtain or print a summary of the results. 
The group testing functions \kbd{predict} (i.e., \code{\link{predict.gt}}) and \kbd{residuals} (i.e., \code{\link{residuals.gt}}) can be used to extract various useful features of the value returned by \kbd{gtreg.halving}. 
An object of class \kbd{"gt"} is a list containing at least the following components: 

 \item{ coefficients }{ a named vector of coefficients }
 \item{ hessian }{ estimated Hessian matrix of the negative log likelihood function, serves as an estimate of the information matrix }
 \item{ residuals }{ the response residuals, difference of the observed group responses and the fitted group responses. }
 \item{ fitted.values }{ the fitted mean values of group responses. }
 \item{ deviance }{ the deviance between the fitted model and the saturated model. }
 \item{ aic }{ Akaike's An Information Criterion, minus twice the maximized log-likelihood plus twice the number of coefficients }
 \item{ counts }{ the number of iterations in the EM algorithm.}
 \item{ z }{ the vector of group responses. }
 \item{ call }{ the matched call. }
 \item{ formula }{ the formula supplied. }
 \item{ terms }{ the terms object used. }
 \item{ method }{ the method \kbd{"Xie"} used to fit the model. }
 \item{ link }{ the link function used in the model. }
}
\references{

\emph{Xie, M. (2001)}, Regression analysis of group testing samples, \emph{Statistics in Medicine, 20, 1957-1969.}

}

\author{
Boan Zhang
}


\seealso{

\code{\link{summary.gt}}, \code{\link{predict.gt}} and \code{\link{residuals.gt}} for \kbd{gt} methods.

}
\examples{

## --- Continuing the Example from  '?sim.halving':

set.seed(46)
gt.data <- sim.halving(par = c(-6, .1), gshape = 17, gscale = 1.4,
               sample.size = 5000, group.size = 5,
               sens = 0.95, spec = 0.95)
fit1 <- gtreg.halving(formula = gres ~ x, data = gt.data, groupn = groupn,
        start = c(-6, .1), subg = subgroup, retest = retest,
        sens = .95, spec = .95, trace = TRUE)
summary(fit1)

}

\keyword{models}
\keyword{regression}
\concept{group testing}
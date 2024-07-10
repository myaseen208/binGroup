\name{gtreg.mp}
\Rdversion{1.1}
\alias{gtreg.mp}
\alias{EM.mp}

\title{
Fitting Group Testing Models in Matrix Pooling Setting
}

\description{
\kbd{gtreg.mp} is a function to fit the group testing regression model in the matrix pooling setting specified through a symbolic description of the linear predictor and descriptions of the group testing setting.
}

\usage{
gtreg.mp(formula, data, coln, rown, arrayn, retest = NULL,
 sens = 1, spec = 1,  linkf = c("logit", "probit", "cloglog"),
 sens.ind = NULL, spec.ind = NULL,  start = NULL,
 control = gt.control(...), ...)

EM.mp(col.resp, row.resp, X, coln, rown, sqn, ret, sens, spec,
 linkf, sens.ind, spec.ind, start = NULL, control = gt.control())

}

\arguments{
  \item{formula}{
an object of class \kbd{"formula"} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under 'Details'.
}
  \item{data}{
an optional data frame, list or environment (or object coercible by \kbd{as.data.frame} to a data frame) containing the variables in the model. If not found in \kbd{data}, the variables are taken from \kbd{environment(formula)}, typically the environment from which \kbd{gtreg.mp} is called.
}
  \item{coln}{
a vector, list or data frame that specifies column group number for each sample
}
  \item{rown}{
a vector, list or data frame that specifies row group number for each sample
}
  \item{arrayn}{
a vector, list or data frame that specifies array number for each sample
}
  \item{retest}{
a vector, list or data frame of individual retest results. A 0 denotes negative and 1 denotes positive. A \code{NA} denotes that no retest is performed for that individual. Default value is \code{NULL} for no retests. 
}
  \item{sens}{
sensitivity of the group tests, set to be 1 by default.
}
  \item{spec}{
specificity of the group tests, set to be 1 by default.
}
  \item{sens.ind}{
sensitivity of the individual retests, set to be equal to \kbd{sens} if not specified otherwise.
}
  \item{spec.ind}{
specificity of the individual retests, set to be equal to \kbd{spec} if not specified otherwise.
}

  \item{linkf}{
a character string specifying one of the three link functions for a binomial model: \kbd{"logit"} (default) or \kbd{"probit"} or \kbd{"cloglog"}.
}

\item{start}{
starting values for the parameters in the linear predictor.
}

\item{control}{a list of parameters for controlling the fitting process. See the documentation for \code{\link{gt.control}} for details.}

  \item{col.resp}{
For \kbd{EM.mp}: vector of group responses of column pools for all samples. 0 denotes negative and 1 denotes positive.
}
  \item{row.resp}{
For \kbd{EM.mp}: vector of group responses of row pools for all samples. 0 denotes negative and 1 denotes positive.
}
  \item{X}{
For \kbd{EM.mp}: the design matrix of the covariates.
}

\item{sqn}{For \kbd{EM.mp}: a vector that specifies array number}

\item{ret}{For \kbd{EM.mp}: a vector containing individual retest results}

\item{\dots}{
arguments to be passed to \code{\link{gt.control}}: see argument \kbd{control}}

}

\details{
With matrix pooling, individual samples are placed in a matrix-like grid where samples are pooled within each row and within each column. This leads to two kinds of group responses: row and column group responses. Thus, a typical predictor has the form \kbd{cbind(col.resp, row.resp)} \kbd{~ covariates} where \kbd{col.resp} is the (numeric) column group response vector and \kbd{row.resp} is the (numeric) row group response vector.
 The covariates term is a series of terms which specifies a linear predictor for individual responses. Note that it is actually the unobserved individual responses, not the observed group responses, which are modeled by the covariates. In \kbd{col.resp} and \kbd{row.resp}, a 0 denotes a negative response and 1 denotes a positive response, where the probability of an individual positive response is being modeled directly. A terms specification of the form \kbd{first + second} indicates all the terms in \kbd{first} together with all the terms in \kbd{second} with duplicates removed. The terms in the formula will be re-ordered so that main effects come first, followed by the interactions, all second-order, all third-order and so on; to avoid this pass a \kbd{terms} object as the formula.

A specification of the form \kbd{first:second} indicates the set of terms obtained by taking the interactions of all terms in \kbd{first} with all terms in \kbd{second}. The specification \kbd{first*second} indicates the cross of \kbd{first} and \kbd{second}. This is the same as \kbd{first + second + first:second}. 

\kbd{EM.mp} is the workhorse function. It applies Xie`s EM algorithm to the likelihood function written in terms of the unobserved individual responses. In each E step, the Gibbs sampling technique is used to estimate the conditional probabilities. Because of the large number of Gibbs samples needed to achieve convergence, the model fitting process could be quite slow, especially when multiple positive rows and columns are observed. In this case, we can either increase the Gibbs sample size to help achieve convergence or loosen the convergence criteria by increasing \kbd{tol} at the expense of perhaps poorer estimates. If follow-up retests are performed, the retest results going into the model will help achieve convergence faster with the same Gibbs sample size and convergence criteria. In each M step, we use \kbd{glm.fit} to update the parameter estimates

}
\value{

\kbd{gtreg.mp} returns an object of class \kbd{"gt.mp"} which inherits from the class \kbd{"gt"}. See later in this section. 
The function \kbd{summary} (i.e., \kbd{summary.gt.mp}) can be used to obtain or print a summary of the results. 
The group testing function \kbd{predict} (i.e., \kbd{predict.gt}) can be used to make predictions on \kbd{"gt.mp"} objects. 
An object of class \kbd{"gt.mp"} is a list containing at least the following components: 

 \item{ coefficients }{ a named vector of coefficients. }
 \item{ hessian }{ estimated Hessian matrix of the negative log likelihood function, serves as an estimate of the information matrix }
 \item{ counts }{ the number of iterations performed in the EM algorithm. }
 \item{ Gibbs.sample.size }{ the number of Gibbs samples generated in each E step. }
 \item{ call }{ the matched call. }
 \item{ formula }{ the formula supplied. }
 \item{ terms }{ the terms object used. }
 \item{ link }{ the link function used in the model. }

}
\references{
\emph{Xie, M. (2001)}, Regression analysis of group testing samples, \emph{Statistics in Medicine, 20, 1957-1969.}
}
\author{
Boan Zhang
}

\seealso{
\code{\link{summary.gt.mp}} and \code{\link{predict.gt}} for \kbd{gt.mp} methods.
\code{\link{gtreg}} for the group testing regression model in the simple pooling setting.

}
\examples{
## --- Continuing the Example from  '?sim.mp':
# 5*6 and 4*5 matrix
set.seed(9128)
sa1a<-sim.mp(par=c(-7,0.1), n.row=c(5,4), n.col=c(6,5),
 sens=0.95, spec=0.95)
sa1<-sa1a$dframe


\dontrun{
fit1 <- gtreg.mp(formula = cbind(col.resp, row.resp) ~ x, data = sa1, 
                 coln = coln, rown = rown, arrayn = arrayn, 
                 sens = 0.95, spec = 0.95, tol = 0.005, n.gibbs = 2000, trace = TRUE)
fit1
summary(fit1)

}

## Here is an example of how long this fitting process may take. For the 
## following simulated data, it takes a computer with 2.2GHZ processor and 
## 3GB RAM about 6 minutes to achieve convergence.
set.seed(9012)
sa2a<-sim.mp(par=c(-7,0.1), n.row=c(10,10,10,10), n.col=c(10,10,10,10), 
             sens=0.95, spec=0.95)
sa2<-sa2a$dframe

\dontrun{
fit2 <- gtreg.mp(formula = cbind(col.resp, row.resp) ~ x, data = sa2, 
                 coln = coln, rown = rown, arrayn = arrayn, retest = retest,
                 sens = 0.95, spec = 0.95, start = c(-7, 0.1), tol = 0.005)

fit2
summary(fit2)

}

}

\keyword{regression}
\concept{matrix pooling}
\concept{group testing}

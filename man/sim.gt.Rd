\name{sim.gt}
\Rdversion{1.1}
\alias{sim.gt}

\title{
Simulation Function for Group Testing Data
}
\description{
Simulates data in group testing form ready to be fit by \kbd{gtreg}.
}
\usage{
sim.gt(x = NULL, gshape = 20, gscale = 2, par,
     linkf = c("logit", "probit", "cloglog"),
     sample.size, group.size, sens = 1, spec = 1,
     sens.ind = NULL, spec.ind = NULL)

}

\arguments{
\item{x}{
a matrix of user-submitted covariates to simulate the data with, defaults to \kbd{NULL} in which case a gamma distribution is used to generate the covariates automatically
}

  \item{gshape}{
shape parameter of gamma distribution, must be non-negative, set to be 20 by default
}
  \item{gscale}{
scale parameter of gamma distribution, must be strictly positive, set to be 2 by default
}
  \item{par}{
the true coefficients in the linear predictor 
}
  \item{sample.size}{
sample size of simulated data 
}

\item{linkf}{
a character string specifying one of the three link functions to be used: \kbd{"logit"} (default) or \kbd{"probit"} or \kbd{"cloglog"}
}

  \item{group.size}{
group size in pooling individual samples
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

}
\details{
\kbd{sim.gt} generates group testing data in simple pooling form. The covariates are either specified by the \kbd{x} argument or they are generated from a gamma distribution with a given \kbd{gshape} and \kbd{gscale}. The individual probabilities are calculated from the covariates, the coefficients given in \kbd{par}, and the link function specified through \kbd{linkf}. The true binary individual responses are then simulated from the individual probabilities. The true group responses are found from the individual responses within the groups (i.e., if at least one response is positive, the group is positive; otherwise, the group response is negative). Finally, the observed group responses are simulated using the given \kbd{sens} and \kbd{spec}. Individual retests are simulated from \kbd{sens.ind} and \kbd{spec.ind} for samples in observed positive groups. Note that with a given \kbd{group.size}, the last group may have less individuals. 

}
\value{
\kbd{sim.gt} returns a data frame with the following columns:
\item{gres}{the group response}
\item{x}{the covariate}
\item{groupn}{the group number}
\item{ind}{the actual individual response}
\item{retest}{the results of individual retests}
}

\author{
Boan Zhang
}

\seealso{
\code{\link{gtreg}}, \code{\link{gtreg.mp}}
}
\examples{

set.seed(46)
gt.data <- sim.gt(par = c(-12, 0.2), sample.size = 700, group.size = 5)

x1 <- sort(runif(100, 0, 30))
x2 <- rgamma(100, shape=17, scale=1.5)
gt.data <- sim.gt(x = cbind(x1,x2), par=c(-14, 0.2, 0.3),
 group.size = 4, sens = 0.98, spec = 0.98)

}

\keyword{datagen}
\concept{group testing}
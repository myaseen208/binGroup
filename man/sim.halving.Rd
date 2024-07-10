\name{sim.halving}
\Rdversion{1.1}
\alias{sim.halving}

\title{
Simulation Function for Group Testing Data for the Halving Protocol
}
\description{
Simulates data in group testing form ready to be fit by \kbd{gtreg.halving}.
}
\usage{
sim.halving(x = NULL, gshape = 20, gscale = 2, par,
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
\kbd{sim.halving} generates group testing data for the halving protocol. The covariates are either specified by the \kbd{x} argument or they are generated from a gamma distribution with a given \kbd{gshape} and \kbd{gscale}. The individual probabilities are calculated from the covariates, the coefficients given in \kbd{par}, and the link function specified through \kbd{linkf}. The true binary individual responses are then simulated from the individual probabilities. The group, subgroup, and individual retests are simulated using the given \kbd{sens} and \kbd{spec} under the halving protocol. 

}
\value{
\kbd{sim.halving} returns a data frame with the following columns:
\item{gres}{the group response}
\item{x}{the covariate}
\item{groupn}{the group number}
\item{ind}{the actual individual response}
\item{retest}{the results of individual retests}
\item{subgroup}{the subgroup number}
}

\author{
Boan Zhang
}

\seealso{
\code{\link{gtreg.halving}}
}
\examples{

set.seed(46)
gt.data <- sim.halving(par = c(-6, .1), gshape = 17, gscale = 1.4,
               sample.size = 5000, group.size = 5,
               sens = 0.95, spec = 0.95)

}

\keyword{datagen}
\concept{group testing}
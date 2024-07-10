\name{sim.mp}
\Rdversion{1.1}
\alias{sim.mp}

\title{
Simulation Function for Group Testing Data with Matrix Pooling Design
}
\description{
Simulates data in group testing form ready to be fit by \kbd{gtreg.mp}.
}
\usage{
sim.mp(x = NULL, gshape = 20, gscale = 2, par,
    linkf = c("logit", "probit", "cloglog"),
    n.row, n.col, sens = 1, spec = 1,
    sens.ind = NULL, spec.ind = NULL)

}

\arguments{
  \item{x}{
a matrix of user-submitted covariates to simulate the data with, defaults to \kbd{NULL} in which case a gamma distribution is used to generate the covariates automatically.
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
\item{linkf}{a character string specifying one of the three link functions to be used: \kbd{"logit"} (default) or \kbd{"probit"} or \kbd{"cloglog"}
}

  \item{n.row}{
a vector that specifies the number of rows in each matrix, a scalar if only one matrix is simulated
}
  \item{n.col}{
a vector that specifies the number of columns in each matrix, a scalar if only one matrix is simulated
}
  \item{sens}{
sensitivity of the group tests, set to be 1 by default.
}
  \item{spec}{
specificity of the group tests, set to be 1 by default.
}
  \item{sens.ind}{
sensitivity of the individual retests, set to be equal to sens if not specified otherwise.
}
  \item{spec.ind}{
specificity of the individual retests, set to be equal to spec if not specified otherwise.
}
}

\details{

\kbd{sim.mp} generates group testing data in matrix pooling form. The covariates are either specified by the \kbd{x} argument or they are generated from a gamma distribution with a given \kbd{gshape} and \kbd{gscale}. The individual probabilities are calculated from the covariates, the coefficients given in \kbd{par} and the link function specified through \kbd{linkf}. The true binary individual responses are then simulated from the individual probabilities. The individuals are organized into (by column) one or more matrices specified by \kbd{n.row} and \kbd{n.col}, and the true group responses are found (i.e., if at least one response is positive, the group is positive; otherwise, the group response is negative). The observed row and column group responses are then simulated using the given \kbd{sens} and \kbd{spec} values. Individual retests are simulated from \kbd{sens.ind} and \kbd{spec.ind} for individuals that lie on the intersection of an observed positive row and an observed positive column. In the case where no column (row) tests positive in a matrix, all the individuals in any observed positive rows (columns) will be assigned a simulated retest result. If no column or row is observed positive, \kbd{NULL} is returned. 

}

\value{

\kbd{sim.mp} returns a list with the components
\kbd{dframe}: the data frame that is actually to be fit,
\kbd{ind}: the true individual responses presented in matrices and
\kbd{prob}: the individual probabilities.

\kbd{dframe} is a data frame with columns
\item{col.resp}{the column group response}
\item{row.resp}{the row group response}
\item{x}{the covariate}
\item{arrayn}{the array number}
\item{coln}{the column group number}
\item{rown}{the row group number}
\item{retest}{the results of individual retests}

}

\author{
Boan Zhang
}

\seealso{
 \code{\link{gtreg.mp}} for the corresponding function to fit the model.
}

\examples{

# 5*6 and 4*5 matrix
set.seed(9128)
sa1a<-sim.mp(par=c(-7,0.1), n.row=c(5,4), n.col=c(6,5),
 sens=0.95, spec=0.95)
sa1<-sa1a$dframe


}

\keyword{datagen}
\concept{group testing}
\concept{matrix pooling}

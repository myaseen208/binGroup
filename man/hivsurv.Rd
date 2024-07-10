\name{hivsurv}
\Rdversion{1.1}
\alias{hivsurv}
\docType{data}

\title{Data from an HIV surveillance project}

\description{
The \kbd{hivsurv} data set comes from an HIV surveillance project discussed in Verstraeten et al. (1998) and Vansteelandt et al. (2000). The purpose of the study was to estimate the HIV prevalence among pregnant Kenyan women in four rural locations of the country, using both individual and group testing responses. Blood tests were administered to each participating woman, and 4 covariates were obtained on each woman. Because the original group responses are unavailable, individuals are artificially put into groups of 5 here to form group responses. Only the 428 complete observations are given here.
}

\usage{data(hivsurv)}

\format{
  A data frame with 428 observations on the following 8 variables.
  \describe{
    \item{\code{DATE}}{the date when each sample was collected}
    \item{\code{PAR.}}{parity (number of children)}
    \item{\code{AGE}}{age (in years)}
    \item{\code{MA.ST.}}{marital status (1: single; 2: married (polygamous); 3: married (monogamous); 4: divorced; 5: widow)}
    \item{\code{EDUC.}}{highest attained education level (1: no schooling; 2: primary school; 3: secondary school; 4: higher)}
    \item{\code{HIV}}{individual response of HIV diagnosis (0: negative; 1: positive)}
    \item{\code{gnum}}{the group number that designates individuals into groups}
    \item{\code{groupres}}{the group response calculated from artificially formed groups}
  }
}

\source{
\emph{Verstraeten, T., Farah, B., Duchateau, L., Matu, R. (1998),} Pooling sera to reduce the cost of HIV surveillance: a feasibility study in a rural Kenyan district, \emph{ Tropical Medicine & International Health, 3, 747-750.}

\emph{Vansteelandt, S., Goetghebeur, E., and Verstraeten, T. (2000),} Regression models for disease prevalence with diagnostic tests on pools of serum samples, \emph{ Biometrics, 56, 1126-1133.}
}

\examples{
data(hivsurv)

str(hivsurv)

}
\keyword{datasets}


# Start  MasterPool.Array.Measures() functions
###################################################################
#    Brianna Hitt - 05-01-17
#    Purpose: calculates descriptive measures for array testing with master pooling
#      calls: f, beta0, beta1, gamma0, gamma1, beta0.y, beta1.y, sum.beta.gamma,
#             mu.y, and nu.y - functions given in the equations for array testing
#             with master pooling in "Comparison of Group Testing Algorithms for Case 
#             Identification in the Presence of Test Error" by Kim et al. (2007)
#      inputs: results = an object containing results from Array.Measures()
#              n = size of a row/column in the square array
#              pmat = matrix of individual probabilities
#              Se = sensitivity of the diagnostic test
#              Sp = specificity of the diagnostic test
#      outputs: list of the expected number of tests (ET), and measures of testing accuracy,
#               including PSe, PSp, PPPV, and PNPV

#' @title Operating characteristics for array testing with master pooling
#' 
#' @description Calculate the expected number of tests and accuracy measures
#' for each individual using array testing with master pooling.
#' 
#' @param results an object containing results (expected number of tests and 
#' accuracy measures) from \code{\link{Array.Measures}}.
#' @param n size of a row/column in the square array.
#' @param pmat matrix of individual risk probabilities.
#' @inheritParams OTC
#' 
#' @details This function assumes that the array is square (i.e., the row
#' and column size are equal) and utilizes the equations from Kim et al. 
#' (2007) for square array testing with master pooling. This function 
#' calculates the operating characteristics for array testing with master 
#' pooling. Operating characteristics calculated are expected number of tests, 
#' pooling sensitivity, pooling specificity, pooling positive predictive value, 
#' and pooling negative predictive value for each individual. 
#'
#' @return A list containing:
#' \item{ET}{the expected number of tests for the array.}
#' \item{PSe}{a matrix containing each individual's pooling sensitivity, 
#' corresponding to the input matrix of individual probabilities.}
#' \item{PSp}{a matrix containing each individual's pooling specificity, 
#' corresponding to the input matrix of individual probabilities.}
#' \item{PPV}{a matrix containing each individual's pooling positive predictive
#' value, corresponding to the input matrix of individual probabilities.}
#' \item{NPV}{a matrix containing each individual's pooling negative predictive
#' value, corresponding to the input matrix of individual probabilities.}
#' 
#' @section Note: This function returns the pooling positive and negative
#' predictive values for all individuals in the array even though these
#' measures are diagnostic specific; i.e., PPV (NPV) should only be considered
#' for those individuals who have tested positive (negative).
#' 
#' @author Brianna D. Hitt
#' 
#' @references 
#' \insertRef{Kim2007}{binGroup}
#' 
#' @seealso
#' \code{\link{Array.Measures}} for calculating operating 
#' characteristics under array testing without master pooling, 
#' \code{\link{hierarchical.desc2}} for three-stage hierarchical and 
#' non-informative two-stage hierarchical testing, and 
#' \code{\link{inf.dorf.measures}} for informative two-stage hierarchical 
#' testing. 
#' 
#' @family Operating characteristic functions
#' 
#' @examples 
#' # Calculate the operating characteristics for 
#' #   non-informative array testing with master
#' #   pooling, with a 6x6 array and an overall 
#' #   disease risk of p = 0.10.
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' p.mat <- matrix(data=0.10, ncol=6, nrow=6)
#' results <- Array.Measures(p=p.mat, se=0.90, sp=0.90)
#' MasterPool.Array.Measures(results=results, n=36, 
#' pmat=p.mat, Se=0.90, Sp=0.90)

MasterPool.Array.Measures <- function(results, n, pmat, Se, Sp){
  
  # extract the measures from the results for Array.Measures()
  ET.A2 <- results$T
  PSe.A2 <- results$PSe
  PSp.A2 <- results$PSp
  PPV.A2 <- results$PPV
  NPV.A2 <- results$NPV
  
  p <- pmat[1,1]
  
  # calculate the measures for array testing with master pooling
  ET.A2M <- 1/(n^2) + (1 - Se - Sp)*((1 - p)^(n^2))*(2/n + (1 - Sp)^2 + 2*(1 - Sp)*Sp^n) + Se*ET.A2
  PSe.A2M <- Se^4 + 2*(Se^3)*(1 - Se)*(1 - f(n=n, p=pmat, Se=Se, Sp=Sp))^(n-1)
  PSp.A2M <- 1 - ((1 - Sp)*mu.y(n=n, p=pmat, Se=Se, Sp=Sp) + 2*(1 - Sp)*nu.y(n=n, p=pmat, Se=Se, Sp=Sp))
  PPV.A2M <- (pmat*PSe.A2M)/((1-pmat)*(1-PSp.A2M) + pmat*PSe.A2M)
  NPV.A2M <- ((1-pmat)*PSp.A2M)/(pmat*(1-PSe.A2M) + (1-pmat)*PSp.A2M)
  
  list("ET"=ET.A2M, "PSe"=PSe.A2M, "PSp"=PSp.A2M, "PPV"=PPV.A2M, "NPV"=NPV.A2M)
}

##################################################################
# Support Functions needed for MasterPool.Array.Measures()       #
##################################################################
f <- function(n, p, Se, Sp){
  (1 - Sp)*(1 - p)^n + Se*(1 - (1 - p)^n)
}

beta0 <- function(n, c, p){
  choose(n=n, k=c)*((1 - p)^(n^2 - n*c + c))*(1 - (1 - p)^(n-1))^c
}

beta1 <- function(n, c, p){
  choose(n=n, k=c)*((1 - p)^(n^2 - n*c))*(1 - (1 - p)^n)^c - beta0(n=n, c=c, p=p)
}

gamma0 <- function(n, c, Se, Sp){
  (1 - Sp)*((1 - Se)^c)*(Sp^(n-c))
}

gamma1 <- function(n, c, Se, Sp){
  Se*(Sp^(n-c))*(1 - Se)^c
}

beta0.y <- function(n, c, p){
  beta0(n=n, c=c, p=p)/(1 - p)
}

beta1.y <- function(n, c, p){
  choose(n=(n-1), k=c)*((1 - (1 - p)^n)^c)*((1 - p)^(n^2 - n*c - 1)) + choose(n=(n-1), k=(c-1))*((1 - (1 - p)^n)^(c-1))*((1 - p)^(n^2 - n*c))*(1 - (1 - p)^(n-1)) - beta0.y(n=n, c=c, p=p)
}

sum.beta.gamma <- function(n,p,Se,Sp){
  sum <- 0
  for(c in 1:n){
    sum <- sum + beta0.y(n=n, c=c, p=p)*gamma0(n=n, c=c, Se=Se, Sp=Sp) + beta1.y(n=n, c=c, p=p)*gamma1(n=n, c=c, Se=Se, Sp=Sp)
  }
  return(sum)
}

mu.y <- function(n, p, Se, Sp){
  f(n=(n^2 - 2*n + 1), p=p, Se=Se, Sp=Sp)*((1 - Sp)*(1 - p)^(n-1))^2 + (Se^2)*(1 - (1 - p)^(n-1))*((1 - Sp)*(1 - p)^(n-1) + f(n=(n - 1), p=p, Se=Se, Sp=Sp))
}

nu.y <- function(n, p, Se, Sp){
  (1 - Sp)*beta0.y(n=n, c=0, p=p)*gamma0(n=n, c=0, Se=Se, Sp=Sp) + Se*sum.beta.gamma(n=n, p=p, Se=Se, Sp=Sp)
}



###################################################################
# Start  NI.A2M() functions
###################################################################

#' @title Find the optimal testing configuration for non-informative 
#' array testing with master pooling
#'
#' @description Find the optimal testing configuration for
#' non-informative array testing with master pooling and
#' calculate the associated operating characteristics.
#'
#' @param p the probability of disease, which can be specified as an overall
#' probability of disease or a homogeneous vector of individual probabilities.
#' @param group.sz a single group size (representing the row/column size)
#' for which to calculate the operating characteristics, or a range of group
#' (row/column) sizes over which to find the OTC.
#' @inheritParams OTC
#'
#' @details This function finds the OTC and computes the associated 
#' operating characteristics for non-informative array testing with 
#' master pooling. Array testing with master pooling involves testing 
#' all specimens in the array together in one group before any row or 
#' column groups are formed. This function uses only square arrays, 
#' which is the way array-based group testing is carried out
#' in most real-world applications. Operating characteristics calculated 
#' are expected number of tests, pooling sensitivity, pooling specificity, 
#' pooling positive predictive value, and pooling negative predictive value 
#' for the algorithm. See Hitt et al. (2018) at
#' \url{http://chrisbilder.com/grouptesting} or Kim et al. (2007) 
#' for additional details on the implementation of non-informative array
#' testing with master pooling.
#' 
#' The value(s) specified by \kbd{group.sz} represent the row/column size, which 
#' is used for the second stage of testing after the entire array is tested 
#' together in one group. If a single value is provided for \kbd{group.sz}, operating
#' characteristics will be calculated and no optimization will be performed.
#' If a range of group sizes is specified, the OTC will be found over all 
#' group sizes.
#' 
#' The displayed pooling sensitivity, pooling specificity, pooling positive 
#' predictive value, and pooling negative predictive value are weighted 
#' averages of the corresponding individual accuracy measures for all 
#' individuals within the initial group for a hierarchical algorithm.
#' Expressions for these averages are provided in the Supplementary 
#' Material for Hitt et al. (2018). These expressions are based on accuracy 
#' definitions given by Altman and Bland (1994a, 1994b).
#'
#' @return A list containing:
#' \item{prob}{the probability of disease, as specified by the user.}
#' \item{Se}{the sensitivity of the diagnostic test.}
#' \item{Sp}{the specificity of the diagnostic test.}
#' \item{opt.ET, opt.MAR, opt.GR}{a list for each objective function specified 
#' by the user, containing:
#' \describe{
#' \item{OTC}{a list specifying elements of the optimal testing configuration, 
#' which include:
#' \describe{
#' \item{Array.dim}{the row/column size for the second stage of testing.}
#' \item{Array.sz}{the overall array size used for the first stage of testing.}}}
#' \item{p.mat}{the matrix of individual probabilities.}
#' \item{ET}{the expected testing expenditure for the OTC.}
#' \item{value}{the value of the objective function per individual.}
#' \item{PSe}{the overall pooling sensitivity for the algorithm.
#' Further details are given under 'Details'.}
#' \item{PSp}{the overall pooling specificity for the algorithm.
#' Further details are given under 'Details'.}
#' \item{PPPV}{the overall pooling positive predictive value for the algorithm.
#' Further details are given under 'Details'.}
#' \item{PNPV}{the overall pooling negative predictive value for the algorithm.
#' Further details are given under 'Details'.}}}
#'
#' @section Note: This function returns the pooling positive and negative
#' predictive values for all individuals in the array even though these
#' measures are diagnostic specific; i.e., PPV (NPV) should only be considered
#' for those individuals who have tested positive (negative).
#' 
#' @author Brianna D. Hitt
#'
#' @references
#' \insertRef{Altman1994a}{binGroup}
#' 
#' \insertRef{Altman1994b}{binGroup}
#' 
#' \insertRef{Graff1972}{binGroup}
#' 
#' \insertRef{Hitt2018}{binGroup}
#' 
#' \insertRef{Kim2007}{binGroup}
#' 
#' \insertRef{Malinovsky2016}{binGroup}
#'
#' @seealso
#' \code{\link{NI.Array}} for non-informative array testing without master 
#' pooling, \code{\link{Inf.Array}} for informative array testing without 
#' master pooling, and \code{\link{OTC}} for finding the optimal testing 
#' configuration for a number of standard group testing algorithms.
#'
#' \url{http://chrisbilder.com/grouptesting}
#'
#' @family OTC functions
#'
#' @examples
#' # Find the OTC for non-informative array testing with 
#' #   master pooling over a range of group (row/column) sizes.
#' # This example takes approximately 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' NI.A2M(p=0.04, Se=0.95, Sp=0.95, group.sz=3:10,
#' obj.fn=c("ET", "MAR"))
#'
#' # Calculate the operating characteristics for a specified 
#' #   group (row/column) size for non-informative array 
#' #   testing with master pooling.
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' NI.A2M(p=rep(0.01, 64), Se=0.90, Sp=0.90, group.sz=8,
#' obj.fn=c("ET", "MAR", "GR"),
#' weights=matrix(data=c(1,1,10,10), 
#' nrow=2, ncol=2, byrow=TRUE))

# Brianna Hitt - 05-01-17
# Updated: Brianna Hitt - 06-20-18
NI.A2M <- function(p, Se, Sp, group.sz, obj.fn, weights=NULL){
  
  start.time<-proc.time()
  
  set.of.I <- group.sz
  
  save.it <- matrix(data=NA, nrow=length(set.of.I), ncol=18)
  count <- 1
  
  for(I in set.of.I){
    N <- I^2
    
    # build a matrix of probabilities
    p.mat <- matrix(data=p[1], nrow=I, ncol=I)
    
    # calculate the measures for array testing, with and without master pooling
    # based on equations provided by Kim et al. (2007) where n represents the row/column size and n^2 represents the array size
    save.info.A2 <- Array.Measures(p=p.mat, se=Se, sp=Sp)
    save.info.A2M <- MasterPool.Array.Measures(results=save.info.A2, n=I, pmat=p.mat, Se=Se, Sp=Sp)
    
    # extract accuracy measures for each individual, with and without master pooling
    ET <- save.info.A2M$ET
    PSe.mat <- save.info.A2M$PSe
    PSp.mat <- save.info.A2M$PSp
    if("MAR" %in% obj.fn){
      MAR <- MAR.func(ET=ET, p.vec=p.mat, PSe.vec=PSe.mat, PSp.vec=PSp.mat)
    } else{MAR <- NA}
    
    # calculate overall accuracy measures
    PSe <- sum(p.mat*PSe.mat)/sum(p.mat)
    PSp <- sum((1-p.mat)*(PSp.mat))/sum(1-p.mat)
    PPPV <- sum(p.mat*PSe.mat)/sum(p.mat*PSe.mat + (1-p.mat)*(1-PSp.mat))
    PNPV <- sum((1-p.mat)*PSp.mat)/sum((1-p.mat)*PSp.mat + p.mat*(1-PSe.mat))
    
    # for each row in the matrix of weights, calculate the GR function
    if(is.null(dim(weights))){
      GR1 <- NA
      GR2 <- NA
      GR3 <- NA
      GR4 <- NA
      GR5 <- NA
      GR6 <- NA
    } else{
      GR1 <- GR.func(ET=ET, p.vec=p.mat, PSe.vec=PSe.mat, PSp.vec=PSp.mat, D1=weights[1,1], D2=weights[1,2])
      if(dim(weights)[1]>=2){
        GR2 <- GR.func(ET=ET, p.vec=p.mat, PSe.vec=PSe.mat, PSp.vec=PSp.mat, D1=weights[2,1], D2=weights[2,2])
      } else{GR2 <- NA}
      if(dim(weights)[1]>=3){
        GR3 <- GR.func(ET=ET, p.vec=p.mat, PSe.vec=PSe.mat, PSp.vec=PSp.mat, D1=weights[3,1], D2=weights[3,2])
      } else{GR3 <- NA}
      if(dim(weights)[1]>=4){
        GR4 <- GR.func(ET=ET, p.vec=p.mat, PSe.vec=PSe.mat, PSp.vec=PSp.mat, D1=weights[4,1], D2=weights[4,2])
      } else{GR4 <- NA}
      if(dim(weights)[1]>=5){
        GR5 <- GR.func(ET=ET, p.vec=p.mat, PSe.vec=PSe.mat, PSp.vec=PSp.mat, D1=weights[5,1], D2=weights[5,2])
      } else{GR5 <- NA}
      if(dim(weights)[1]>=6){
        GR6 <- GR.func(ET=ET, p.vec=p.mat, PSe.vec=PSe.mat, PSp.vec=PSp.mat, D1=weights[6,1], D2=weights[6,2])
      } else{GR6 <- NA}
    }
    
    save.it[count,] <- c(p[1], Se, Sp, I, N, ET, ET/N, MAR, GR1/N, GR2/N, GR3/N, GR4/N, GR5/N, GR6/N, PSe, PSp, PPPV, PNPV)
    cat("Row/Column Size=", I, ", Array Size=", N, "\n", sep="")
    count <- count + 1
  }
  
  # find the optimal testing configuration, over all array sizes considered
  result.ET <- save.it[save.it[,7]==min(save.it[,7]), c(1:6,7,15:ncol(save.it))]
  result.MAR <- save.it[save.it[,8]==min(save.it[,8]), c(1:6,8,15:ncol(save.it))]
  result.GR1 <- save.it[save.it[,9]==min(save.it[,9]), c(1:6,9,15:ncol(save.it))]
  result.GR2 <- save.it[save.it[,10]==min(save.it[,10]), c(1:6,10,15:ncol(save.it))]
  result.GR3 <- save.it[save.it[,11]==min(save.it[,11]), c(1:6,11,15:ncol(save.it))]
  result.GR4 <- save.it[save.it[,12]==min(save.it[,12]), c(1:6,12,15:ncol(save.it))]
  result.GR5 <- save.it[save.it[,13]==min(save.it[,13]), c(1:6,13,15:ncol(save.it))]
  result.GR6 <- save.it[save.it[,14]==min(save.it[,14]), c(1:6,14,15:ncol(save.it))]
  
  p.mat.ET <- matrix(data=result.ET[1], nrow=result.ET[4], ncol=result.ET[4])
  if("MAR" %in% obj.fn){
    p.mat.MAR <- matrix(data=result.MAR[1], nrow=result.MAR[4], ncol=result.MAR[4])
  } else{p.mat.MAR <- NA}
  if(is.null(dim(weights))){
    p.mat.GR1 <- NA
    p.mat.GR2 <- NA
    p.mat.GR3 <- NA
    p.mat.GR4 <- NA
    p.mat.GR5 <- NA
    p.mat.GR6 <- NA
  } else{
    p.mat.GR1 <- matrix(data=result.GR1[1], nrow=result.GR1[4], ncol=result.GR1[4])
    if(dim(weights)[1]>=2){
      p.mat.GR2 <- matrix(data=result.GR2[1], nrow=result.GR2[4], ncol=result.GR2[4])
    } else{p.mat.GR2 <- NA}
    if(dim(weights)[1]>=3){
      p.mat.GR3 <- matrix(data=result.GR3[1], nrow=result.GR3[4], ncol=result.GR3[4])
    } else{p.mat.GR3 <- NA}
    if(dim(weights)[1]>=4){
      p.mat.GR4 <- matrix(data=result.GR4[1], nrow=result.GR4[4], ncol=result.GR4[4])
    } else{p.mat.GR4 <- NA}
    if(dim(weights)[1]>=5){
      p.mat.GR5 <- matrix(data=result.GR5[1], nrow=result.GR5[4], ncol=result.GR5[4])
    } else{p.mat.GR5 <- NA}
    if(dim(weights)[1]>=6){
      p.mat.GR6 <- matrix(data=result.GR6[1], nrow=result.GR6[4], ncol=result.GR6[4])
    } else{p.mat.GR6 <- NA}
  }
  
  # create a list of results for each objective function
  opt.ET <- list("OTC"=list("Array.dim"=result.ET[4], "Array.sz"=result.ET[5]), "p.mat"=p.mat.ET, "ET"=result.ET[6], "value"=result.ET[7], "PSe"=result.ET[8], "PSp"=result.ET[9], "PPPV"=result.ET[10], "PNPV"=result.ET[11])
  opt.MAR <- list("OTC"=list("Array.dim"=result.MAR[4], "Array.sz"=result.MAR[5]), "p.mat"=p.mat.MAR, "ET"=result.MAR[6], "value"=result.MAR[7], "PSe"=result.MAR[8], "PSp"=result.MAR[9], "PPPV"=result.MAR[10], "PNPV"=result.MAR[11])
  opt.GR1 <- list("OTC"=list("Array.dim"=result.GR1[4], "Array.sz"=result.GR1[5]), "p.mat"=p.mat.GR1, "ET"=result.GR1[6], "value"=result.GR1[7], "PSe"=result.GR1[8], "PSp"=result.GR1[9], "PPPV"=result.GR1[10], "PNPV"=result.GR1[11])
  opt.GR2 <- list("OTC"=list("Array.dim"=result.GR2[4], "Array.sz"=result.GR2[5]), "p.mat"=p.mat.GR2, "ET"=result.GR2[6], "value"=result.GR2[7], "PSe"=result.GR2[8], "PSp"=result.GR2[9], "PPPV"=result.GR2[10], "PNPV"=result.GR2[11])
  opt.GR3 <- list("OTC"=list("Array.dim"=result.GR3[4], "Array.sz"=result.GR3[5]), "p.mat"=p.mat.GR3, "ET"=result.GR3[6], "value"=result.GR3[7], "PSe"=result.GR3[8], "PSp"=result.GR3[9], "PPPV"=result.GR3[10], "PNPV"=result.GR3[11])
  opt.GR4 <- list("OTC"=list("Array.dim"=result.GR4[4], "Array.sz"=result.GR4[5]), "p.mat"=p.mat.GR4, "ET"=result.GR4[6], "value"=result.GR4[7], "PSe"=result.GR4[8], "PSp"=result.GR4[9], "PPPV"=result.GR4[10], "PNPV"=result.GR4[11])
  opt.GR5 <- list("OTC"=list("Array.dim"=result.GR5[4], "Array.sz"=result.GR5[5]), "p.mat"=p.mat.GR5, "ET"=result.GR5[6], "value"=result.GR5[7], "PSe"=result.GR5[8], "PSp"=result.GR5[9], "PPPV"=result.GR5[10], "PNPV"=result.GR5[11])
  opt.GR6 <- list("OTC"=list("Array.dim"=result.GR6[4], "Array.sz"=result.GR6[5]), "p.mat"=p.mat.GR6, "ET"=result.GR6[6], "value"=result.GR6[7], "PSe"=result.GR6[8], "PSp"=result.GR6[9], "PPPV"=result.GR6[10], "PNPV"=result.GR6[11])
  
  # create a list of results, including all objective functions
  opt.all <- list("opt.ET"=opt.ET, "opt.MAR"=opt.MAR, "opt.GR1"=opt.GR1, "opt.GR2"=opt.GR2, 
                  "opt.GR3"=opt.GR3, "opt.GR4"=opt.GR4, "opt.GR5"=opt.GR5, "opt.GR6"=opt.GR6)
  # remove any objective functions not requested by the user
  opt.req <- Filter(function(x) !is.na(x$ET), opt.all)
  time.it(start.time)
  c("prob"=list(p), "Se"=Se, "Sp"=Sp, opt.req)
  
}

###################################################################

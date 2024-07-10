#' @title Find the optimal testing configuration
#' 
#' @description Find the optimal testing configuration (OTC) for standard group testing
#' algorithms and calculate the associated operating characteristics.
#'
#' @param algorithm character string defining the group testing algorithm to be used.
#' Non-informative testing options include two-stage hierarchical ("\kbd{D2}"),
#' three-stage hierarchical ("\kbd{D3}"), square array testing without master 
#' pooling ("\kbd{A2}"), and square array testing without master pooling ("\kbd{A2M}"). 
#' Informative testing options include two-stage hierarchical ("\kbd{ID2}"),
#' three-stage hierarchical ("\kbd{ID3}"), and square array testing without 
#' master pooling ("\kbd{IA2}").
#' @param p overall probability of disease that will be used to generate a
#' vector/matrix of individual probabilities. For non-informative algorithms, a 
#' homogeneous set of probabilities will be used. For informative algorithms, the 
#' \code{\link{p.vec.func}} function will be used to generate a heterogeneous 
#' set of probabilities. Either \kbd{p} or \kbd{probabilities} should be specified, 
#' but not both.
#' @param probabilities a vector of individual probabilities, which is homogeneous 
#' for non-informative testing algorithms and heterogeneous for informative 
#' testing algorithms. Either  \kbd{p} or \kbd{probabilities} should be specified, 
#' but not both.
#' @param Se the sensitivity of the diagnostic test.
#' @param Sp the specificity of the diagnostic test.
#' @param group.sz a single group size or range of group sizes for which to 
#' calculate operating characteristics and/or find the OTC. The details of group 
#' size specification are given under 'Details'. 
#' @param obj.fn a list of objective functions which are minimized to find the
#' OTC. The expected number of tests per individual, "\kbd{ET}", will always 
#' be calculated. Additional options include "\kbd{MAR}" 
#' (the expected number of tests divided by the expected number of correct 
#' classifications, described in Malinovsky et al. (2016)), and "\kbd{GR}" 
#' (a linear combination of the expected number of tests, the number of 
#' misclassified negatives, and the number of misclassified positives, 
#' described in Graff & Roeloffs (1972)). See Hitt et al. (2018) at
#' \url{http://chrisbilder.com/grouptesting} for additional details.
#' @param weights a matrix of up to six sets of weights for the GR function. 
#' Each set of weights is specified by a row of the matrix.
#' @param alpha a shape parameter for the beta distribution that specifies the degree of
#' heterogeneity for the generated probability vector (for informative testing only).
#'
#' @details This function finds the OTC and computes the
#' associated operating characteristics for standard group testing algorithms,
#' as described in Hitt et al. (2018) at 
#' \url{http://chrisbilder.com/grouptesting}.
#'
#' Available algorithms include two- and three-stage hierarchical testing and
#' array testing with and without master pooling. Both non-informative and informative
#' group testing settings are allowed for each algorithm, except informative 
#' array testing with master pooling is unavailable because this method has not 
#' appeared in the group testing literature. Operating characteristics calculated are
#' expected number of tests, pooling sensitivity, pooling specificity, pooling
#' positive predictive value, and pooling negative predictive value for each individual.
#'
#' The value(s) specified by \kbd{group.sz} represent the initial (stage 1) 
#' group size for three-stage hierarchical testing and non-informative 
#' two-stage hierarchical testing. For informative two-stage hierarchical testing, 
#' the \kbd{group.sz} specified represents the block size used in the pool-specific
#' optimal Dorfman (PSOD) method, where the initial group (block) is not
#' tested. For more details on informative two-stage hierarchical testing 
#' implemented via the PSOD method, see Hitt et al. (2018) and McMahan et al. (2012a).
#' For array testing without master pooling, the \kbd{group.sz} specified
#' represents the row/column size for initial (stage 1) testing. For array testing 
#' with master pooling, the \kbd{group.sz} specified represents the row/column size 
#' for stage 2 testing. The group size for initial testing is overall array size, 
#' given by the square of the row/column size.
#' 
#' If a single value is provided for \kbd{group.sz} with array testing or
#' non-informative two-stage hierarchical testing, operating
#' characteristics will be calculated and no optimization will be performed.
#' If a single value is provided for \kbd{group.sz} with three-stage hierarchical or 
#' informative two-stage hierarchical, the OTC will be  
#' found over all possible configurations. If a range of group sizes is specified, 
#' the OTC will be found over all group sizes.
#' 
#' The displayed pooling sensitivity, pooling specificity, pooling positive 
#' predictive value, and pooling negative predictive value are weighted 
#' averages of the corresponding individual accuracy measures for all 
#' individuals within the initial group for a hierarchical algorithm, or 
#' within the entire array for an array-based algorithm.
#' Expressions for these averages are provided in the Supplementary 
#' Material for Hitt et al. (2018). These expressions are based on accuracy 
#' definitions given by Altman and Bland (1994a, 1994b).
#' 
#' @return A list containing:
#' \item{prob}{the probability of disease, as specified by the user.}
#' \item{alpha}{level of heterogeneity for the generated probability vector
#' (for informative testing only).}
#' \item{Se}{the sensitivity of the diagnostic test.}
#' \item{Sp}{the specificity of the diagnostic test.}
#' \item{opt.ET, opt.MAR, opt.GR}{a list for each objective function specified
#' by the user, containing:
#' \describe{
#' \item{OTC}{a list specifying elements of the optimal testing configuration, 
#' which may include:
#' \describe{
#' \item{Stage1}{pool size for the first stage of hierarchical testing, if applicable.}
#' \item{Stage2}{pool sizes for the second stage of hierarchical testing, if applicable.}
#' \item{Block.sz}{the block size/initial group size for informative Dorfman testing,
#' which is not tested.}
#' \item{pool.szs}{pool sizes for the first stage of testing for informative Dorfman
#' testing.}
#' \item{Array.dim}{the row/column size for array testing.}
#' \item{Array.sz}{the overall array size for array testing (the square of the row/column size).}}}
#' \item{p.vec}{the sorted vector of individual probabilities, if applicable.}
#' \item{p.mat}{the sorted matrix of individual probabilities in gradient arrangement,
#' if applicable.}
#' \item{ET}{the expected testing expenditure for the OTC.}
#' \item{value}{the value of the objective function per individual.}
#' \item{PSe}{the overall pooling sensitivity for the algorithm. Further details 
#' are given under 'Details'.}
#' \item{PSp}{the overall pooling specificity for the algorithm. Further details 
#' are given under 'Details'.}
#' \item{PPPV}{the overall pooling positive predictive value for the algorithm. 
#' Further details are given under 'Details'.}
#' \item{PNPV}{the overall pooling negative predictive value for the algorithm. 
#' Further details are given under 'Details'.}}}
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
#' \insertRef{Malinovsky2016}{binGroup}
#' 
#' \insertRef{McMahan2012a}{binGroup}
#' 
#' \insertRef{McMahan2012b}{binGroup}
#' 
#'
#' @seealso
#' \code{\link{NI.Dorf}} for non-informative two-stage (Dorfman) testing, \code{\link{Inf.Dorf}} for
#' informative two-stage (Dorfman) testing, \code{\link{NI.D3}} for non-informative three-stage
#' hierarchical testing, \code{\link{Inf.D3}} for informative three-stage hierarchical testing,
#' \code{\link{NI.Array}} for non-informative array testing, \code{\link{Inf.Array}} for informative
#' array testing, and \code{\link{NI.A2M}} for non-informative array testing with master pooling.
#'
#' \url{http://chrisbilder.com/grouptesting}
#'
#' @family OTC functions
#' 
#' @examples
#' # Find the OTC for non-informative
#' #   two-stage hierarchical (Dorfman) testing
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' OTC(algorithm="D2", p=0.05, Se=0.99, Sp=0.99, group.sz=2:100,
#' obj.fn=c("ET", "MAR"))
#'
#' # Find the OTC for informative
#' #   two-stage hierarchical (Dorfman) testing, implemented
#' #   via the pool-specific optimal Dorfman (PSOD) method
#' #   described in McMahan et al. (2012a), where the greedy
#' #   algorithm proposed for PSOD is replaced by considering
#' #   all possible testing configurations.
#' # A vector of individual probabilities is generated using
#' #   the expected value of order statistics from a beta 
#' #   distribution with p = 0.01 and a heterogeneity level 
#' #   of alpha = 0.5. Depending on the specified probability, 
#' #   alpha level, and overall group size, simulation may 
#' #   be necessary in order to generate the vector of individual
#' #   probabilities. This is done using p.vec.func() and 
#' #   requires the user to set a seed in order to reproduce 
#' #   results.
#' # This example takes approximately 2.5 minutes to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' \dontrun{
#' set.seed(52613)
#' OTC(algorithm="ID2", p=0.01, Se=0.95, Sp=0.95, group.sz=50,
#' obj.fn=c("ET", "MAR", "GR"),
#' weights=matrix(data=c(1, 1, 10, 10, 0.5, 0.5), 
#' nrow=3, ncol=2, byrow=TRUE), alpha=0.5)}
#'
#' # Find the OTC over all possible
#' #   testing configurations for a specified group size for
#' #   non-informative three-stage hierarchical testing
#' # This example takes approximately 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' OTC(algorithm="D3", p=0.001, Se=0.95, Sp=0.95, group.sz=18,
#' obj.fn=c("ET", "MAR", "GR"),
#' weights=matrix(data=c(1, 1), nrow=1, ncol=2, byrow=TRUE))
#'
#' # Find the OTC for non-informative
#' #   three-stage hierarchical testing
#' # This example takes approximately 20 seconds to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' \dontrun{
#' OTC(algorithm="D3", p=0.06, Se=0.90, Sp=0.90, 
#' group.sz=3:30, obj.fn=c("ET", "MAR", "GR"),
#' weights=matrix(data=c(1, 1, 10, 10, 100, 100), 
#' nrow=3, ncol=2, byrow=TRUE))}
#'
#' # Find the OTC over all possible configurations 
#' #   for a specified group size, given a 
#' #   heterogeneous vector of probabilities.
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' OTC(algorithm="ID3", probabilities=c(0.012, 0.014, 0.011, 
#' 0.012, 0.010, 0.015), Se=0.99, Sp=0.99, group.sz=6, 
#' obj.fn=c("ET","MAR","GR"), weights=matrix(data=c(1, 1), 
#' nrow=1, ncol=2, byrow=TRUE), alpha=0.5)
#'
#' # Calculate the operating characteristics for a specified array size
#' #   for non-informative array testing without master pooling
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' OTC(algorithm="A2", p=0.005, Se=0.95, Sp=0.95, group.sz=8, 
#' obj.fn=c("ET", "MAR"))
#'
#' # Find the OTC for informative array testing without 
#' #   master pooling
#' # A vector of individual probabilities is generated using
#' #   the expected value of order statistics from a beta 
#' #   distribution with p = 0.03 and a heterogeneity level 
#' #   of alpha = 2. The probabilities are then arranged in 
#' #   a matrix using the gradient method described in 
#' #   McMahan et al. (2012b). Depending on the specified 
#' #   probability, alpha level, and overall group size, 
#' #   simulation may be necessary in order to generate the 
#' #   vector of individual probabilities. This is done using 
#' #   p.vec.func() and requires the user to set a 
#' #   seed in order to reproduce results.
#' # This example takes approximately 30 seconds to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' \dontrun{
#' set.seed(1002)
#' OTC(algorithm="IA2", p=0.03, Se=0.95, Sp=0.95, 
#' group.sz=3:20, obj.fn=c("ET", "MAR", "GR"),
#' weights=matrix(data=c(1, 1, 10, 10, 100, 100), 
#' nrow=3, ncol=2, byrow=TRUE), alpha=2)}
#'
#' # Find the OTC for non-informative array testing 
#' #   with master pooling
#' # This example takes approximately 20 seconds to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' \dontrun{
#' OTC(algorithm="A2M", p=0.02, Se=0.90, Sp=0.90, 
#' group.sz=3:20, obj.fn=c("ET", "MAR", "GR"),
#' weights=matrix(data=c(1, 1, 10, 10, 0.5, 0.5, 2, 2, 
#' 100, 100, 10, 100), nrow=6, ncol=2, byrow=TRUE))}

OTC <- function(algorithm, p=NULL, probabilities=NULL, Se=0.99, Sp=0.99, group.sz, obj.fn=c("ET","MAR"), weights=NULL, alpha=2){
  
  ## make sure that all necessary information is included in the correct format
  if(is.null(p) & is.null(probabilities)){
    stop("Please specify an overall probability of disease using the 'p' argument, \n or specify a vector of individual probabilities using the 'probabilities' argument.")
  } else if(!is.null(p) & !is.null(probabilities)){
    stop("You have specified both an overall probability of disease AND a \n vector of individual probabilities. Please specify only one option.")
  } else{
    if(!is.null(p)){
      if(length(p)==1){
        cat("You have specified an overall probability of disease. \n A probability vector will be generated based on the algorithm specified.\n")
      } else{
        stop("You have specified a probability vector instead of an overall probability of disease.\n Please specify an overall probability of disease, and the probability vector will be \n generated based on the algorithm specified for each group size included in the range.\n")
      }
    }
    if(!is.null(probabilities)){
      if(length(group.sz)==1){
        cat("You have specified a vector containing individual probabilities of disease.\n")
        if((algorithm %in% c("D2", "D3", "ID2", "ID3")) & length(probabilities)!=group.sz){
          stop("The vector of individual probabilities is not the correct length. Please make sure\n that the length of the probability vector is the same as the specified group size.\n")
        } else if((algorithm %in% c("A2", "A2M", "IA2")) & length(probabilities)!=group.sz^2){
          stop("The vector of individual probabilities is not the correct length. Please make sure that the\n length of the probability vector is the same as the overall array size (the square of the specified group size).\n")
        }
        if((algorithm %in% c("D2", "D3", "A2", "A2M")) & all.equal(probabilities, rep(probabilities[1],length(probabilities)))!=TRUE){
          stop("You have specified a heterogeneous probability vector for a non-informative\n algorithm. Please specify a homogeneous probability vector using the 'probabilities'\n argument or specify an overall probability of disease using the 'p' argument.\n")
        }
      } else if(length(group.sz)>1){
        stop("You have specified a probability vector along with a range \n of group sizes. Please specify a single group size.\n")
      }
    }
  }
  
  if(length(group.sz)==1){
    if(algorithm %in% c("D3", "ID2", "ID3")){
      cat("A single group size was provided. The optimal testing configuration will be found \n over all possible testing configurations for the specified group size.\n")
    } else{
      cat("A single group size was provided. No optimization will be performed.\n")
    }
  }
  
  if(is.null(obj.fn)){
    stop("Please specify one or more objective functions for which to find the optimal testing configuration.\n")
  }
  
  if("GR" %in% obj.fn){
    if(is.null(weights)){
      stop("No weights have been specified. The GR function will not be calculated.\n")
    } else if(dim(weights)[2]!=2){
      stop("Please check the dimension of the weights matrix. \n Each row should specify a set of weights, D1 and D2.\n")
    }
  }
  
  # call function for non-informative two-stage hierarchical (Dorfman) testing
  if(algorithm == "D2"){
    if(min(group.sz)<2){
      stop("Please specify a minimum group size of at least 2.\n")
    }
    if(!is.null(alpha)){
      cat("Non-informative Testing - alpha will be ignored. \n")
    }
    cat("Algorithm: Non-informative two-stage hierarchical (Dorfman) testing \n")
    if(!is.null(p)){
      results <- NI.Dorf(p=p, Se=Se, Sp=Sp, group.sz=group.sz, obj.fn=obj.fn, weights=weights)
    } else if(!is.null(probabilities)){
      results <- NI.Dorf(p=probabilities, Se=Se, Sp=Sp, group.sz=group.sz, obj.fn=obj.fn, weights=weights)
    }
  }
  
  # call function for non-informative three-stage hierarchical testing
  if(algorithm == "D3"){
    if(min(group.sz)<3){
      stop("Please specify a minimum group size of at least 3.\n")
    }
    if(max(group.sz)>50){
      message("NOTE: You have specified a maximum group size of 50 or larger.\n This function may take a VERY long time to run.\n Press 'ESC' if you wish to cancel the submitted statements.\n")
    }
    if(!is.null(alpha)){
      cat("Non-informative Testing - alpha will be ignored. \n")
    }
    cat("Algorithm: Non-informative three-stage hierarchical testing \n")
    if(!is.null(p)){
      results <- NI.D3(p=p, Se=Se, Sp=Sp, group.sz=group.sz, obj.fn=obj.fn, weights=weights)
    } else if(!is.null(probabilities)){
      results <- NI.D3(p=probabilities, Se=Se, Sp=Sp, group.sz=group.sz, obj.fn=obj.fn, weights=weights)
    }
  }
  
  # call function for non-informative square array testing without master pooling
  if(algorithm == "A2"){
    if(min(group.sz)<3){
      stop("Please specify a minimum group size of at least 3.\n")
    }
    if(max(group.sz)>50){
      message("NOTE: You have specified a maximum group size of 50 or larger.\n This function may take a VERY long time to run.\n Press 'ESC' if you wish to cancel the submitted statements.\n")
    }
    if(!is.null(alpha)){
      cat("Non-informative Testing - alpha will be ignored. \n")
    }
    cat("Algorithm: Non-informative square array testing without master pooling \n")
    if(!is.null(p)){
      results <- NI.Array(p=p, Se=Se, Sp=Sp, group.sz=group.sz, obj.fn=obj.fn, weights=weights)
    } else if(!is.null(probabilities)){
      results <- NI.Array(p=probabilities, Se=Se, Sp=Sp, group.sz=group.sz, obj.fn=obj.fn, weights=weights)
    }
  }
  
  # call function for non-informative square array testing with master pooling
  if(algorithm == "A2M"){
    if(min(group.sz)<3){
      stop("Please specify a minimum group size of at least 3.\n")
    }
    if(max(group.sz)>50){
      message("NOTE: You have specified a maximum group size of 50 or larger.\n This function may take a VERY long time to run.\n Press 'ESC' if you wish to cancel the submitted statements.\n")
    }
    if(!is.null(alpha)){
      cat("Non-informative Testing - alpha will be ignored. \n")
    }
    cat("Algorithm: Non-informative square array testing with master pooling \n")
    if(!is.null(p)){
      results <- NI.A2M(p=p, Se=Se, Sp=Sp, group.sz=group.sz, obj.fn=obj.fn, weights=weights)
    } else if(!is.null(probabilities)){
      results <- NI.A2M(p=probabilities, Se=Se, Sp=Sp, group.sz=group.sz, obj.fn=obj.fn, weights=weights)
    }
  }
  
  # call function for informative two-stage hierarchical (Dorfman) testing
  if(algorithm == "ID2"){
    if(min(group.sz)<2){
      stop("Please specify a minimum group size of at least 2.\n")
    }
    if(max(group.sz)>=50){
      message("NOTE: You have specified a maximum group size of 50 or larger.\n This function may take a VERY long time to run.\n Press 'ESC' if you wish to cancel the submitted statements.\n")
    }
    cat("Algorithm: Informative Dorfman testing \n")
    if(!is.null(p)){
      results <- Inf.Dorf(p=p, Se=Se, Sp=Sp, group.sz=group.sz, obj.fn=obj.fn, weights=weights, alpha=alpha)
    } else if(!is.null(probabilities)){
      results <- Inf.Dorf(p=probabilities, Se=Se, Sp=Sp, group.sz=group.sz, obj.fn=obj.fn, weights=weights, alpha=alpha)
    }
  }
  
  # call function for informative three-stage hierarchical testing
  if(algorithm == "ID3"){
    if(min(group.sz)<3){
      stop("Please specify a minimum group size of at least 3.\n")
    }
    if(max(group.sz)>=50){
      message("NOTE: You have specified a maximum group size of 50 or larger.\n This function may take a VERY long time to run.\n Press 'ESC' if you wish to cancel the submitted statements.\n")
    }
    cat("Algorithm: Informative three-stage hierarchical testing \n")
    if(!is.null(p)){
      results <- Inf.D3(p=p, Se=Se, Sp=Sp, group.sz=group.sz, obj.fn=obj.fn, weights=weights, alpha=alpha)
    } else if(!is.null(probabilities)){
      results <- Inf.D3(p=probabilities, Se=Se, Sp=Sp, group.sz=group.sz, obj.fn=obj.fn, weights=weights, alpha=alpha)
    }
  }
  
  # call function for informative square array testing without master pooling
  if(algorithm == "IA2"){
    if(min(group.sz)<3){
      stop("Please specify a minimum group size of at least 3.\n")
    }
    if(max(group.sz)>=50){
      message("NOTE: You have specified a maximum group size of 50 or larger.\n This function may take a VERY long time to run.\n Press 'ESC' if you wish to cancel the submitted statements.\n")
    }
    cat("Algorithm: Informative square array testing without master pooling \n")
    if(!is.null(p)){
      results <- Inf.Array(p=p, Se=Se, Sp=Sp, group.sz=group.sz, obj.fn=obj.fn, weights=weights, alpha=alpha)
    } else if(!is.null(probabilities)){
      results <- Inf.Array(p=probabilities, Se=Se, Sp=Sp, group.sz=group.sz, obj.fn=obj.fn, weights=weights, alpha=alpha)
    }
  }
  
  results
}



# Supporting functions for OTC and the associated calls
####################################################################

#' @title Generate a vector of probabilities for informative group 
#' testing algorithms.
#' 
#' @description Generate a vector of individual risk probabilities
#' using an overall probability of disease (i.e., the expected 
#' value of order statistics from a beta distribution)
#' for use with informative group testing algorithms.
#' 
#' @param p overall probability of disease that will be used to 
#' generate a vector of individual risk probabilities.
#' @param alpha a shape parameter for the beta distribution that
#' specifies the degree of heterogeneity for the generated
#' probability vector.
#' @param grp.sz the number of total individuals for which to 
#' generate risk probabilities.
#' 
#' @details This function uses Michael Black's \code{\link{beta.dist}}
#' function to generate a vector of individual risk probabilities,
#' ordered from least to greatest. Depending on the specified 
#' probability, alpha level, and overall group size, simulation 
#' may be necessary in order to generate the vector of individual 
#' probabilities. For this reason, the user should set a seed in
#' order to reproduce results. The \kbd{p.vec.func} function
#' augments the \code{\link{beta.dist}} function by checking whether 
#' simulation is needed before attempting to generate the vector
#' of individual risk probabilities. See Black et al. (2015)
#' for additional details on Michael Black's \code{\link{beta.dist}}
#' function.
#' 
#' @return A vector of individual risk probabilities.
#' 
#' @author Brianna D. Hitt
#' 
#' @references 
#' \insertRef{Black2015}{binGroup}
#' 
#' @seealso 
#' \code{\link{Informative.array.prob}} for arranging a vector
#' of individual risk probabilities in a matrix for informative
#' array testing without master pooling and \code{\link{beta.dist}}
#' for the function on which \kbd{p.vec.func} is based on.
#' 
#' \url{http://chrisbilder.com/grouptesting}
#' 
#' @family Individual risk probability functions
#' 
#' @examples
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' set.seed(8791)
#' p.vec.func(p=0.03, alpha=0.5, grp.sz=100)
#' 
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' set.seed(52613)
#' p.vec.func(p=0.005, alpha=2, grp.sz=40)

p.vec.func <- function(p, alpha, grp.sz){
  if(is.na(p)){
    NA
  } else{
    p.try <- suppressWarnings(try(beta.dist(p=p, alpha=alpha, grp.sz=grp.sz), silent=TRUE)) 
    if(class(p.try)=="try-error"){
      beta.dist(p=p, alpha=alpha, grp.sz=grp.sz, simul=TRUE)
    } else{
      beta.dist(p=p, alpha=alpha, grp.sz=grp.sz)
    } 
  }
}

###################################################################




# Start MAR.func() function
###################################################################
#    Brianna Hitt - 4-17-17
#    Purpose: calculates MAR objective function, from Malinovsky, Albert & Roy (2015)
#      inputs: ET - expected number of tests
#              p.vec - vector of individual probabilities
#              PSe.vec - vector of individual pooling sensitivities
#              PSp.vec - vector of individual pooling specificities
#      note: The MAR objective function divides ET, the expected number of tests, by EC,
#            the expected number of correct classifications, and should be minimized.
#            Note: Malinovsky, Albert, & Roy (2015) maximized the reciprocal, E(C)/E(T).

MAR.func <- function(ET, p.vec, PSe.vec, PSp.vec){
  EC <- sum(PSe.vec*p.vec + PSp.vec*(1-p.vec))
  ET/EC
}
###################################################################




# Start GR.func() function
###################################################################
#    Brianna Hitt - 4-17-17
#    Purpose: calculates GR objective function, from Graff & Roeloffs (1972)
#               M = E(T) + D_1*(# of misclassified negatives) + D_2*(# of misclassified positives)
#      inputs: ET - expected number of tests
#              p.vec - vector of individual probabilities
#              PSe.vec - vector of individual pooling sensitivities
#              PSp.vec - vector of individual pooling specificities
#              D1, D2 - weights/costs for misclassification
#      note: this function specifies equal weights of 1 by default

GR.func <- function(ET, p.vec, PSe.vec, PSp.vec, D1=1, D2=1){
  ET + D1*sum((1-PSp.vec)*(1-p.vec)) + D2*sum((1-PSe.vec)*p.vec)
}
###################################################################




# Start time.it() function
###################################################################
#    Brianna Hitt - 5-13-17
#    Purpose: calculates the time elapsed
#      inputs: x = object containing the start time

time.it <- function(x) {
  end.time<-proc.time()
  save.time<-end.time-x
  cat("\n Number of minutes running:", save.time[3]/60, "\n \n")
  save.time[3]/60
}
###################################################################

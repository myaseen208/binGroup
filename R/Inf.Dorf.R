
# Start  inf.dorf.measures() function
###################################################################
#    Brianna Hitt - 4-18-17
#    Purpose: calculates the testing expenditure and accuracy measures for informative Dorfman testing,
#             given a testing configuration; Informative Dorfman testing is the same as pool-specific
#             optimal Dorfman testing, described by McMahan et al. (2012), except for that it attempts
#             to find the optimal testing configuration by examining all possible configurations rather
#             than using the greedy algorithm proposed by McMahan et al. (2012)
#      calls: characteristics.pool() - calculates the testing expenditure for a given configuration
#             accuracy.dorf() - calculates measures of testing accuracy for informative Dorfman testing
#      inputs: p = probability that an individual is infected, can be an overall probability
#                  of disease or a vector of individual probabilities
#              Se = sensitivity of the diagnostic test
#              Sp = specificity of the diagnostic test
#              N = block size/initial group  size, up to 50
#              pool.sizes = a configuration/set of pool sizes from a matrix of all
#                    possible configurations, generated using the three-stage hierarchical setup
#      outputs: list of the testing expenditure (e), testing variance (v), and
#               measures of testing accuracy (summary)
#      notes: much of the code for this function, including the characteristics.pool() and
#             accuracy.dorf() functions are from Chris McMahan's programs, provided with
#             "Informative Dorfman screening" by McMahan, Tebbs, and Bilder (2012)

#' @title Operating characteristics for informative two-stage 
#' hierarchical (Dorfman) testing
#' 
#' @description Calculate the expected number of tests and accuracy 
#' measures for each individual using informative two-stage 
#' hierarchical (Dorfman) testing, given a vector of individual 
#' probabilities and a testing configuration.
#' 
#' @param prob vector of probabilities corresponding to each individual's 
#' risk of disease.
#' @param N block size/initial group size that is not tested. This
#' is the total number of individuals being tested. The 
#' details of \kbd{block.sz} are given under 'Details'.
#' @param pool.sizes a vector of pool sizes for the first stage 
#' of testing.
#' @inheritParams hierarchical.desc2
#' 
#' @details This function utilizes the equations given by McMahan et al. 
#' (2012) for informative two-stage hierarchical (Dorfman) testing. 
#' It also repurposes functions written by Christopher S. McMahan (see
#' \url{http://chrisbilder.com/grouptesting})
#' for the implementation of informative Dorfman testing and directly 
#' uses functions written for the calculation of the associated operating 
#' characteristics. This function calculates the operating characteristics 
#' for informative two-stage hierarchical (Dorfman) testing. Operating 
#' characteristics calculated are expected number of tests, and pooling 
#' sensitivity, pooling specificity, pooling positive predictive value, 
#' and pooling negative predictive value for each individual.
#' 
#' The specified \kbd{N} represents the block size used in the pool-specific
#' optimal Dorfman (PSOD) method. This is the total number of individuals being
#' tested by the algorithm. This block is not initially tested. Instead, 
#' multiple initial pool sizes within this block are found and tested in 
#' the first stage of testing. The second stage of testing consists of individual
#' retesting. For more information on block size specification, see McMahan et
#' al. (2012).
#' 
#' @return A list containing:
#' \item{e}{the expected number of tests needed to decode all N
#' individuals.}
#' \item{v}{the variance of the total number of tests needed to 
#' decode all N individuals.}
#' \item{summary}{a matrix containing the pool, probability of 
#' disease, pooling sensitivity, pooling specificity, pooling 
#' positive predictive value, and pooling negative predictive 
#' value for each individual. The pool column identifies which 
#' pool each individual is contained in for the first stage of
#' testing.}
#' 
#' @section Note: This function returns the pooling positive and negative
#' predictive values for all individuals even though these measures are 
#' diagnostic specific; i.e., PPPV (PNPV) should only be considered
#' for those individuals who have tested positive (negative).
#' 
#' @author The majority of this function was originally written by 
#' Christopher S. McMahan for McMahan et al. (2012). The function was 
#' obtained from \url{http://chrisbilder.com/grouptesting}. Minor modifications 
#' were made to the function for inclusion in the binGroup package.
#' 
#' @references 
#' \insertRef{McMahan2012a}{binGroup}
#' 
#' @seealso
#' \code{\link{Array.Measures}} for calculating operating characteristics
#' under array testing without master pooling, 
#' \code{\link{MasterPool.Array.Measures}} for non-informative array 
#' testing with master pooling, and \code{\link{inf.dorf.measures}} 
#' for informative two-stage hierarchical testing. See 
#' \code{\link{p.vec.func}} for generating a vector of 
#' individual risk probabilities for informative group testing. 
#' 
#' This function repurposed code from \code{\link{opt.info.dorf}} so that
#' PSOD testing could be implemented using all possible testing configurations
#' instead of a greedy algorithm. This function also used
#' \code{\link{characteristics.pool}} and \code{\link{accuracy.dorf}}
#' to calculate operating characteristics for the optimal set of 
#' pool sizes.
#' 
#' \url{http://chrisbilder.com/grouptesting}
#' 
#' @family Operating characteristic functions
#' @family Informative Dorfman functions
#' 
#' @examples 
#' # Calculate the operating characteristics for 
#' #   informative two-stage hierarchical (Dorfman) testing
#' #   with an overall disease prevalence of E(p(i)) = 0.01,
#' #   where a block size of 50 is split into initial pools
#' #   of 18, 13, 11, and 8 individuals.
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' set.seed(8791)
#' inf.dorf.measures(prob=beta.dist(p=0.01, alpha=2, 
#' grp.sz=50, simul=TRUE), se=0.95, sp=0.95, N=50, 
#' pool.sizes=c(18, 13, 11, 8))

inf.dorf.measures <- function(prob, se, sp, N, pool.sizes){
  
  # Saves original ordering
  ind.order<-order(prob)
  
  # Orders subjects, required under all Informative algorithms
  prob<-sort(prob)
  
  # Determines the number of subjects being screened and sets up vectors for storing summary measures
  N<-length(prob)
  pool.id<-rep(-1000, N)
  PSe<-rep(-100, N)
  PSp<-rep(-100, N)
  PPV<-rep(-100, N)
  NPV<-rep(-100, N)
  
  # pool.sizes is a single configuration/set of pool sizes from the matrix of all
  # possible configurations from the three-stage hierarchical testing setup
  psz <- pool.sizes
  J <- length(psz)
  
  # Finds measures pool by pool
  psz<-c(psz, 0)
  lower<-1
  upper<-psz[1]
  vec.e<-rep(-1000, J)
  vec.v<-rep(-1000, J)
  for(i in 1:J){
    p.pool<-prob[lower:upper]
    pool.id[lower:upper]<-rep(i, length(p.pool))
    
    # calculates the testing expenditure and variance per individual
    res<-characteristics.pool(p=p.pool, se=se, sp=sp)
    vec.e[i]<-res$e
    vec.v[i]<-res$v
    
    # calculates measures of testing accuracy per individual
    res.acc<-accuracy.dorf(p=p.pool, se=se, sp=sp)
    PSe[lower:upper]<-res.acc$PSe
    PSp[lower:upper]<-res.acc$PSp
    PPV[lower:upper]<-res.acc$PPV
    NPV[lower:upper]<-res.acc$NPV
    
    lower<-1+upper
    upper<-upper+psz[i+1]
  }
  
  # Finds total expectation and variation
  res.e<-sum(vec.e)
  res.v<-sum(vec.v)
  
  # Returns all subjects to original ordering, along with their corresponding measures
  prob<-prob[order(ind.order)]
  pool.id<-pool.id[order(ind.order)]
  PSe<-PSe[order(ind.order)]
  PSp<-PSp[order(ind.order)]
  PPV<-PPV[order(ind.order)]
  NPV<-NPV[order(ind.order)]
  
  res.mat<-matrix(c(pool.id, prob, PSe, PSp, PPV, NPV), nrow=N, ncol=6, byrow=FALSE, 
                  dimnames=list(as.character(1:N) , c("pool", "probability", "PSe", "PSp", "PPV", "NPV")))
  
  prob<-prob[ind.order]
  
  list("e"=res.e, "v"=res.v, "summary"=res.mat)
}

###################################################################




# Start  Inf.Dorf() function
###################################################################

#' @title Find the optimal testing configuration for informative 
#' two-stage hierarchical (Dorfman) testing
#'
#' @description Find the optimal testing configuration (OTC) for 
#' informative two-stage hierarchical (Dorfman) testing and 
#' calculate the associated operating characteristics.
#'
#' @param p the probability of disease, which can be specified as an overall 
#' probability of disease, from which a heterogeneous vector of individual 
#' probabilities will be generated, or a heterogeneous vector of individual 
#' probabilities specified by the user.
#' @param group.sz a single block size for which to find the OTC
#' out of all possible configurations, or a range of block sizes over 
#' which to find the OTC.
#' @param alpha a scale parameter for the beta distribution that specifies 
#' the degree of heterogeneity for the generated probability vector. If a 
#' heterogeneous vector of individual probabilities is specified by the user, 
#' \kbd{alpha} can be specified as \kbd{NA} or will be ignored.
#' @inheritParams OTC
#'
#' @details This function finds the OTC and computes the associated operating 
#' characteristics for informative two-stage hierarchical (Dorfman) testing, 
#' implemented via the pool-specific optimal Dorfman (PSOD) method described in 
#' McMahan et al. (2012). This function finds the optimal testing configuration 
#' by considering all possible testing configurations instead of using the greedy 
#' algorithm proposed for PSOD testing. Operating characteristics calculated are
#' expected number of tests, pooling sensitivity, pooling specificity, pooling
#' positive predictive value, and pooling negative predictive value for the algorithm.
#' See Hitt et al. (2018) or McMahan et al. (2012) at 
#' \url{http://chrisbilder.com/grouptesting} for additional details on the 
#' implementation of informative two-stage hierarchical (Dorfman) testing.
#'
#' The value(s) specified by \kbd{group.sz} represent the overall block size 
#' used in the pool-specific optimal Dorfman (PSOD) method, where the overall group 
#' size is not tested. Instead, multiple initial pool sizes within this block are 
#' found and tested in the first stage of testing. The second stage of testing consists 
#' of individual retesting. For more details on informative two-stage hierarchical testing 
#' implemented via the PSOD method, see Hitt et al. (2018) and McMahan et al. (2012).
#' 
#' If a single value is provided for \kbd{group.sz}, the OTC will be  
#' found over all possible testing configurations. If a range of group sizes is 
#' specified, the OTC will be found over all group sizes.
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
#' \item{alpha}{the level of heterogeneity used to generate the vector of individual
#' probabilities.}
#' \item{Se}{the sensitivity of the diagnostic test.}
#' \item{Sp}{the specificity of the diagnostic test.}
#' \item{opt.ET, opt.MAR, opt.GR}{a list for each objective function specified 
#' by the user, containing:
#' \describe{
#' \item{OTC}{a list specifying elements of the optimal testing configuration, 
#' which include:
#' \describe{
#' \item{Block.sz}{the block size/overall group size, which is not tested.}
#' \item{pool.szs}{pool sizes for the first stage of testing.}}}
#' \item{p.vec}{the sorted vector of individual probabilities.}
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
#' @author Brianna D. Hitt
#'
#' @references
#' \insertRef{Altman1994a}{binGroup}
#' 
#' \insertRef{Altman1994b}{binGroup}
#' 
#' \insertRef{Dorfman1943}{binGroup}
#' 
#' \insertRef{Graff1972}{binGroup}
#' 
#' \insertRef{Hitt2018}{binGroup}
#' 
#' \insertRef{Malinovsky2016}{binGroup}
#' 
#' \insertRef{McMahan2012a}{binGroup}
#'
#' @seealso
#' \code{\link{NI.Dorf}} for non-informative two-stage hierarchical (Dorfman) 
#' testing and \code{\link{OTC}} for finding the optimal testing configuration 
#' for a number of standard group testing algorithms.
#'
#' \url{http://chrisbilder.com/grouptesting}
#'
#' @family OTC functions
#'
#' @examples
#' # Find the OTC for informative two-stage hierarchical 
#' #   (Dorfman) testing.
#' # A vector of individual probabilities is generated using
#' #   the expected value of order statistics from a beta 
#' #   distribution with p = 0.01 and a heterogeneity level 
#' #   of alpha = 2. Depending on the specified probability, 
#' #   alpha level, and overall group size, simulation may 
#' #   be necessary in order to generate the vector of individual
#' #   probabilities. This is done using p.vec.func() and 
#' #   requires the user to set a seed in order to reproduce 
#' #   results.
#' # This example takes approximately 20 seconds to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' \dontrun{
#' set.seed(9245)
#' Inf.Dorf(p=0.01, Se=0.95, Sp=0.95, group.sz=3:30, 
#' obj.fn=c("ET", "MAR"), alpha=2)}
#' 
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' set.seed(9245)
#' Inf.Dorf(p=0.01, Se=0.95, Sp=0.95, group.sz=5:10, 
#' obj.fn=c("ET", "MAR"), alpha=2)
#'
#' # Find the OTC for informative two-stage hierarchical 
#' #   (Dorfman) testing, for a specified block size.
#' # This example uses rbeta() to generate random probabilities 
#' #   and requires the user to set a seed in order to reproduce 
#' #   results.
#' # This example takes approximately 2.5 minutes to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' \dontrun{
#' set.seed(8791)
#' Inf.Dorf(p=p.vec.func(p=0.03, alpha=0.5, grp.sz=50), 
#' Se=0.90, Sp=0.90, group.sz=50, obj.fn=c("ET", "MAR", "GR"),
#' weights=matrix(data=c(1,1,10,10), nrow=2, ncol=2, byrow=TRUE),
#' alpha=NA)}

# Brianna Hitt - 4-18-17
# Updated: Brianna Hitt - 06-20-18
Inf.Dorf <- function(p, Se, Sp, group.sz, obj.fn, weights=NULL, alpha=2){
  
  start.time<-proc.time()
  
  set.of.blocks <- group.sz
  
  save.ET <- matrix(data=NA, nrow=length(set.of.blocks), ncol=2*max(set.of.blocks)+10)
  save.MAR <- matrix(data=NA, nrow=length(set.of.blocks), ncol=2*max(set.of.blocks)+10)
  save.GR1 <- matrix(data=NA, nrow=length(set.of.blocks), ncol=2*max(set.of.blocks)+10)
  save.GR2 <- matrix(data=NA, nrow=length(set.of.blocks), ncol=2*max(set.of.blocks)+10)
  save.GR3 <- matrix(data=NA, nrow=length(set.of.blocks), ncol=2*max(set.of.blocks)+10)
  save.GR4 <- matrix(data=NA, nrow=length(set.of.blocks), ncol=2*max(set.of.blocks)+10)
  save.GR5 <- matrix(data=NA, nrow=length(set.of.blocks), ncol=2*max(set.of.blocks)+10)
  save.GR6 <- matrix(data=NA, nrow=length(set.of.blocks), ncol=2*max(set.of.blocks)+10)
  
  count <- 1
  
  for(N in set.of.blocks){
    # build a vector of probabilities for a heterogeneous population
    if(length(p)==1){
      p.vec <- p.vec.func(p=p, alpha=alpha, grp.sz=N)
    } else if(length(p)>1){
      p.vec <- sort(p)
      alpha <- NA
    }
    
    # generate a matrix of all possible configurations/sets of pool sizes
    # the parts() function requires loading the partitions library
    # do not include the first column, which would test the initial group twice
    possible.groups <- parts(n=N)[,-1]
    
    save.it <- matrix(data=NA, nrow=ncol(possible.groups), ncol=2*max(set.of.blocks)+17)
    
    counter <- 1
    for(c in 1:ncol(possible.groups)){
      pool.sizes <- possible.groups[,c]
      # calculate descriptive measures for informative Dorfman testing, given a configuration/set of pool sizes
      save.info <- inf.dorf.measures(prob=p.vec, se=Se, sp=Sp, N=N, pool.sizes=pool.sizes[pool.sizes!=0])
      
      # extract the configuration/pool sizes
      row.names(save.info$summary)=NULL
      pool.sz <- table(save.info$summary[,1])
      row.names(pool.sz)=NULL
      
      # extract accuracy measures for each individual
      ET <- save.info$e
      PSe.vec <- save.info$summary[,3]
      PSp.vec <- save.info$summary[,4]
      if("MAR" %in% obj.fn){
        MAR <- MAR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec)
      } else{MAR <- NA}
      
      # calculate overall accuracy measures
      PSe <- sum(p.vec*PSe.vec)/sum(p.vec)
      PSp <- sum((1-p.vec)*(PSp.vec))/sum(1-p.vec)
      PPPV <- sum(p.vec*PSe.vec)/sum(p.vec*PSe.vec + (1-p.vec)*(1-PSp.vec))
      PNPV <- sum((1-p.vec)*PSp.vec)/sum((1-p.vec)*PSp.vec + p.vec*(1-PSe.vec))
      
      # for each row in the matrix of weights, calculate the GR function
      if(is.null(dim(weights))){
        GR1 <- NA
        GR2 <- NA
        GR3 <- NA
        GR4 <- NA
        GR5 <- NA
        GR6 <- NA
      } else{
        GR1 <- GR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec, D1=weights[1,1], D2=weights[1,2])
        if(dim(weights)[1]>=2){
          GR2 <- GR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec, D1=weights[2,1], D2=weights[2,2])
        } else{GR2 <- NA}
        if(dim(weights)[1]>=3){
          GR3 <- GR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec, D1=weights[3,1], D2=weights[3,2])
        } else{GR3 <- NA}
        if(dim(weights)[1]>=4){
          GR4 <- GR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec, D1=weights[4,1], D2=weights[4,2])
        } else{GR4 <- NA}
        if(dim(weights)[1]>=5){
          GR5 <- GR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec, D1=weights[5,1], D2=weights[5,2])
        } else{GR5 <- NA}
        if(dim(weights)[1]>=6){
          GR6 <- GR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec, D1=weights[6,1], D2=weights[6,2])
        } else{GR6 <- NA}
      }
      
      save.it[counter,] <- c(p.vec, rep(NA, max(0, max(set.of.blocks)-length(p.vec))), alpha, Se, Sp, N, ET, ET/N, MAR, GR1/N, GR2/N, GR3/N, GR4/N, GR5/N, GR6/N, PSe, PSp, PPPV, PNPV, pool.sz, rep(0, max(0, max(set.of.blocks)-length(pool.sz))))
      counter <- counter + 1
      
    }
    
    # find the best configuration for each block size N, out of all possible configurations
    save.ET[count,] <- save.it[save.it[,(max(set.of.blocks)+6)]==min(save.it[,(max(set.of.blocks)+6)]), c(1:(max(set.of.blocks)+5),(max(set.of.blocks)+6),(max(set.of.blocks)+14):ncol(save.it))]
    if(class(try(save.MAR[count,] <- save.it[save.it[,(max(set.of.blocks)+7)]==min(save.it[,(max(set.of.blocks)+7)]), c(1:(max(set.of.blocks)+5),(max(set.of.blocks)+7),(max(set.of.blocks)+14):ncol(save.it))],silent=T))!="try-error"){
      save.MAR[count,] <- save.it[save.it[,(max(set.of.blocks)+7)]==min(save.it[,(max(set.of.blocks)+7)]), c(1:(max(set.of.blocks)+5),(max(set.of.blocks)+7),(max(set.of.blocks)+14):ncol(save.it))]
    }
    if(class(try(save.GR1[count,] <- save.it[save.it[,(max(set.of.blocks)+8)]==min(save.it[,(max(set.of.blocks)+8)]), c(1:(max(set.of.blocks)+5),(max(set.of.blocks)+8),(max(set.of.blocks)+14):ncol(save.it))],silent=T))!="try-error"){
      save.GR1[count,] <- save.it[save.it[,(max(set.of.blocks)+8)]==min(save.it[,(max(set.of.blocks)+8)]), c(1:(max(set.of.blocks)+5),(max(set.of.blocks)+8),(max(set.of.blocks)+14):ncol(save.it))]
    }
    if(class(try(save.GR2[count,] <- save.it[save.it[,(max(set.of.blocks)+9)]==min(save.it[,(max(set.of.blocks)+9)]), c(1:(max(set.of.blocks)+5),(max(set.of.blocks)+9),(max(set.of.blocks)+14):ncol(save.it))],silent=T))!="try-error"){
      save.GR2[count,] <- save.it[save.it[,(max(set.of.blocks)+9)]==min(save.it[,(max(set.of.blocks)+9)]), c(1:(max(set.of.blocks)+5),(max(set.of.blocks)+9),(max(set.of.blocks)+14):ncol(save.it))]
    }
    if(class(try(save.GR3[count,] <- save.it[save.it[,(max(set.of.blocks)+10)]==min(save.it[,(max(set.of.blocks)+10)]), c(1:(max(set.of.blocks)+5),(max(set.of.blocks)+10),(max(set.of.blocks)+14):ncol(save.it))],silent=T))!="try-error"){
      save.GR3[count,] <- save.it[save.it[,(max(set.of.blocks)+10)]==min(save.it[,(max(set.of.blocks)+10)]), c(1:(max(set.of.blocks)+5),(max(set.of.blocks)+10),(max(set.of.blocks)+14):ncol(save.it))]
    }
    if(class(try(save.GR4[count,] <- save.it[save.it[,(max(set.of.blocks)+11)]==min(save.it[,(max(set.of.blocks)+11)]), c(1:(max(set.of.blocks)+5),(max(set.of.blocks)+11),(max(set.of.blocks)+14):ncol(save.it))],silent=T))!="try-error"){
      save.GR4[count,] <- save.it[save.it[,(max(set.of.blocks)+11)]==min(save.it[,(max(set.of.blocks)+11)]), c(1:(max(set.of.blocks)+5),(max(set.of.blocks)+11),(max(set.of.blocks)+14):ncol(save.it))]
    }
    if(class(try( save.GR5[count,] <- save.it[save.it[,(max(set.of.blocks)+12)]==min(save.it[,(max(set.of.blocks)+12)]), c(1:(max(set.of.blocks)+5),(max(set.of.blocks)+12),(max(set.of.blocks)+14):ncol(save.it))],silent=T))!="try-error"){
      save.GR5[count,] <- save.it[save.it[,(max(set.of.blocks)+12)]==min(save.it[,(max(set.of.blocks)+12)]), c(1:(max(set.of.blocks)+5),(max(set.of.blocks)+12),(max(set.of.blocks)+14):ncol(save.it))]
    }
    if(class(try(save.GR6[count,] <- save.it[save.it[,(max(set.of.blocks)+13)]==min(save.it[,(max(set.of.blocks)+13)]), c(1:(max(set.of.blocks)+5),(max(set.of.blocks)+13),(max(set.of.blocks)+14):ncol(save.it))],silent=T))!="try-error"){
      save.GR6[count,] <- save.it[save.it[,(max(set.of.blocks)+13)]==min(save.it[,(max(set.of.blocks)+13)]), c(1:(max(set.of.blocks)+5),(max(set.of.blocks)+13),(max(set.of.blocks)+14):ncol(save.it))]    
    }
    
    cat("Block Size =", N, "\n")
    count <- count + 1
  }
  
  # find the optimal configuration over all block sizes considered
  result.ET <- save.ET[save.ET[,(max(set.of.blocks)+6)]==min(save.ET[,(max(set.of.blocks)+6)]),]
  result.MAR <- save.MAR[save.MAR[,(max(set.of.blocks)+6)]==min(save.MAR[,(max(set.of.blocks)+6)]),]
  result.GR1 <- save.GR1[save.GR1[,(max(set.of.blocks)+6)]==min(save.GR1[,(max(set.of.blocks)+6)]),]
  result.GR2 <- save.GR2[save.GR2[,(max(set.of.blocks)+6)]==min(save.GR2[,(max(set.of.blocks)+6)]),]
  result.GR3 <- save.GR3[save.GR3[,(max(set.of.blocks)+6)]==min(save.GR3[,(max(set.of.blocks)+6)]),]
  result.GR4 <- save.GR4[save.GR4[,(max(set.of.blocks)+6)]==min(save.GR4[,(max(set.of.blocks)+6)]),]
  result.GR5 <- save.GR5[save.GR5[,(max(set.of.blocks)+6)]==min(save.GR5[,(max(set.of.blocks)+6)]),]
  result.GR6 <- save.GR6[save.GR6[,(max(set.of.blocks)+6)]==min(save.GR6[,(max(set.of.blocks)+6)]),]
  
  p.vec.ET <- (result.ET[1:max(set.of.blocks)])[!is.na(result.ET[1:max(set.of.blocks)])]
  if("MAR" %in% obj.fn){
    p.vec.MAR <- (result.MAR[1:max(set.of.blocks)])[!is.na(result.MAR[1:max(set.of.blocks)])]
  } else{p.vec.MAR <- NA}
  if(is.null(dim(weights))){
    p.vec.GR1 <- NA
    p.vec.GR2 <- NA
    p.vec.GR3 <- NA
    p.vec.GR4 <- NA
    p.vec.GR5 <- NA
    p.vec.GR6 <- NA
  } else{
    p.vec.GR1 <- (result.GR1[1:max(set.of.blocks)])[!is.na(result.GR1[1:max(set.of.blocks)])]
    if(dim(weights)[1]>=2){
      p.vec.GR2 <- (result.GR2[1:max(set.of.blocks)])[!is.na(result.GR2[1:max(set.of.blocks)])]
    } else{p.vec.GR2 <- NA}
    if(dim(weights)[1]>=3){
      p.vec.GR3 <- (result.GR3[1:max(set.of.blocks)])[!is.na(result.GR3[1:max(set.of.blocks)])]
    } else{p.vec.GR3 <- NA}
    if(dim(weights)[1]>=4){
      p.vec.GR4 <- (result.GR4[1:max(set.of.blocks)])[!is.na(result.GR4[1:max(set.of.blocks)])]
    } else{p.vec.GR4 <- NA}
    if(dim(weights)[1]>=5){
      p.vec.GR5 <- (result.GR5[1:max(set.of.blocks)])[!is.na(result.GR5[1:max(set.of.blocks)])]
    } else{p.vec.GR5 <- NA}
    if(dim(weights)[1]>=6){
      p.vec.GR6 <- (result.GR6[1:max(set.of.blocks)])[!is.na(result.GR6[1:max(set.of.blocks)])]
    } else{p.vec.GR6 <- NA}
  }
  
  # create a list of results for each objective function
  opt.ET <- list("OTC"=list("Block.sz"=result.ET[(max(set.of.blocks)+4)], "pool.szs"=(result.ET[(max(set.of.blocks)+11):length(result.ET)])[result.ET[(max(set.of.blocks)+11):length(result.ET)]!=0]), "p.vec"=p.vec.ET, 
                 "ET"=result.ET[(max(set.of.blocks)+5)], "value"=result.ET[(max(set.of.blocks)+6)], "PSe"=result.ET[(max(set.of.blocks)+7)], "PSp"=result.ET[(max(set.of.blocks)+8)], "PPPV"=result.ET[(max(set.of.blocks)+9)], "PNPV"=result.ET[(max(set.of.blocks)+10)])
  opt.MAR <- list("OTC"=list("Block.sz"=result.MAR[(max(set.of.blocks)+4)], "pool.szs"=(result.MAR[(max(set.of.blocks)+11):length(result.MAR)])[result.MAR[(max(set.of.blocks)+11):length(result.MAR)]!=0]), "p.vec"=p.vec.MAR,
                  "ET"=result.MAR[(max(set.of.blocks)+5)], "value"=result.MAR[(max(set.of.blocks)+6)], "PSe"=result.MAR[(max(set.of.blocks)+7)], "PSp"=result.MAR[(max(set.of.blocks)+8)], "PPPV"=result.MAR[(max(set.of.blocks)+9)], "PNPV"=result.MAR[(max(set.of.blocks)+10)])  
  opt.GR1 <- list("OTC"=list("Block.sz"=result.GR1[(max(set.of.blocks)+4)], "pool.szs"=(result.GR1[(max(set.of.blocks)+11):length(result.GR1)])[result.GR1[(max(set.of.blocks)+11):length(result.GR1)]!=0]), "p.vec"=p.vec.GR1,
                  "ET"=result.GR1[(max(set.of.blocks)+5)], "value"=result.GR1[(max(set.of.blocks)+6)], "PSe"=result.GR1[(max(set.of.blocks)+7)], "PSp"=result.GR1[(max(set.of.blocks)+8)], "PPPV"=result.GR1[(max(set.of.blocks)+9)], "PNPV"=result.GR1[(max(set.of.blocks)+10)])
  opt.GR2 <- list("OTC"=list("Block.sz"=result.GR2[(max(set.of.blocks)+4)], "pool.szs"=(result.GR2[(max(set.of.blocks)+11):length(result.GR2)])[result.GR2[(max(set.of.blocks)+11):length(result.GR2)]!=0]), "p.vec"=p.vec.GR2,
                  "ET"=result.GR2[(max(set.of.blocks)+5)], "value"=result.GR2[(max(set.of.blocks)+6)], "PSe"=result.GR2[(max(set.of.blocks)+7)], "PSp"=result.GR2[(max(set.of.blocks)+8)], "PPPV"=result.GR2[(max(set.of.blocks)+9)], "PNPV"=result.GR2[(max(set.of.blocks)+10)])  
  opt.GR3 <- list("OTC"=list("Block.sz"=result.GR3[(max(set.of.blocks)+4)], "pool.szs"=(result.GR3[(max(set.of.blocks)+11):length(result.GR3)])[result.GR3[(max(set.of.blocks)+11):length(result.GR3)]!=0]), "p.vec"=p.vec.GR3,
                  "ET"=result.GR3[(max(set.of.blocks)+5)], "value"=result.GR3[(max(set.of.blocks)+6)], "PSe"=result.GR3[(max(set.of.blocks)+7)], "PSp"=result.GR3[(max(set.of.blocks)+8)], "PPPV"=result.GR3[(max(set.of.blocks)+9)], "PNPV"=result.GR3[(max(set.of.blocks)+10)])
  opt.GR4 <- list("OTC"=list("Block.sz"=result.GR4[(max(set.of.blocks)+4)], "pool.szs"=(result.GR4[(max(set.of.blocks)+11):length(result.GR4)])[result.GR4[(max(set.of.blocks)+11):length(result.GR4)]!=0]), "p.vec"=p.vec.GR4,
                  "ET"=result.GR4[(max(set.of.blocks)+5)], "value"=result.GR4[(max(set.of.blocks)+6)], "PSe"=result.GR4[(max(set.of.blocks)+7)], "PSp"=result.GR4[(max(set.of.blocks)+8)], "PPPV"=result.GR4[(max(set.of.blocks)+9)], "PNPV"=result.GR4[(max(set.of.blocks)+10)])
  opt.GR5 <- list("OTC"=list("Block.sz"=result.GR5[(max(set.of.blocks)+4)], "pool.szs"=(result.GR5[(max(set.of.blocks)+11):length(result.GR5)])[result.GR5[(max(set.of.blocks)+11):length(result.GR5)]!=0]), "p.vec"=p.vec.GR5,
                  "ET"=result.GR5[(max(set.of.blocks)+5)], "value"=result.GR5[(max(set.of.blocks)+6)], "PSe"=result.GR5[(max(set.of.blocks)+7)], "PSp"=result.GR5[(max(set.of.blocks)+8)], "PPPV"=result.GR5[(max(set.of.blocks)+9)], "PNPV"=result.GR5[(max(set.of.blocks)+10)])
  opt.GR6 <- list("OTC"=list("Block.sz"=result.GR6[(max(set.of.blocks)+4)], "pool.szs"=(result.GR6[(max(set.of.blocks)+11):length(result.GR6)])[result.GR6[(max(set.of.blocks)+11):length(result.GR6)]!=0]), "p.vec"=p.vec.GR6,
                  "ET"=result.GR6[(max(set.of.blocks)+5)], "value"=result.GR6[(max(set.of.blocks)+6)], "PSe"=result.GR6[(max(set.of.blocks)+7)], "PSp"=result.GR6[(max(set.of.blocks)+8)], "PPPV"=result.GR6[(max(set.of.blocks)+9)], "PNPV"=result.GR6[(max(set.of.blocks)+10)])
  
  # create a list of results, including all objective functions
  opt.all <- list("opt.ET"=opt.ET, "opt.MAR"=opt.MAR, "opt.GR1"=opt.GR1, "opt.GR2"=opt.GR2, 
                  "opt.GR3"=opt.GR3, "opt.GR4"=opt.GR4, "opt.GR5"=opt.GR5, "opt.GR6"=opt.GR6)
  # remove any objective functions not requested by the user
  opt.req <- Filter(function(x) !is.na(x$ET), opt.all)
  time.it(start.time)
  c("prob"=list(p), "alpha"=alpha, "Se"=Se, "Sp"=Sp, opt.req)
  
}

###################################################################


# Start  Inf.D3() function
###################################################################

#' @title Find the optimal testing configuration for informative 
#' three-stage hierarchical testing
#'
#' @description Find the optimal testing configuration (OTC) for 
#' informative three-stage hierarchical testing and calculate the 
#' associated operating charcteristics.
#'
#' @param p the probability of disease, which can be specified as an overall 
#' probability of disease, from which a heterogeneous vector of individual 
#' probabilities will be generated, or a heterogeneous vector of individual 
#' probabilities specified by the user.
#' @param group.sz a single group size over which to find the OTC
#' out of all possible testing configurations, or a range of group sizes
#' over which to find the OTC.
#' @param alpha a scale parameter for the beta distribution that specifies 
#' the degree of heterogeneity for the generated probability vector. If a 
#' heterogeneous vector of individual probabilities is specified by the user, 
#' \kbd{alpha} can be specified as \kbd{NA} or will be ignored.
#' @inheritParams OTC
#'
#' @details This function finds the OTC and computes the associated 
#' operating characteristics for informative three-stage hierarchical testing. 
#' This function finds the optimal testing configuration by considering all 
#' possible testing configurations. Operating characteristics calculated are
#' expected number of tests, pooling sensitivity, pooling specificity, pooling
#' positive predictive value, and pooling negative predictive value for the algorithm.
#' See Hitt et al. (2018) or Black et al. (2015) at 
#' \url{http://chrisbilder.com/grouptesting} for additional details on the 
#' implementation of informative three-stage hierarchical testing.
#' 
#' The value(s) specified by \kbd{group.sz} represent the initial (stage 1) 
#' group size. If a single value is provided for \kbd{group.sz}, the OTC 
#' will be found over all possible testing configurations for that initial 
#' group size. If a range of group sizes is specified, the OTC will be found 
#' over all group sizes.
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
#' \item{Stage1}{pool size for the first stage of testing, i.e. the initial
#' group size.}
#' \item{Stage2}{pool sizes for the second stage of testing.}}}
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
#' \insertRef{Black2015}{binGroup}
#' 
#' \insertRef{Graff1972}{binGroup}
#' 
#' \insertRef{Hitt2018}{binGroup}
#' 
#' \insertRef{Malinovsky2016}{binGroup}
#'
#' @seealso
#' \code{\link{NI.D3}} for non-informative three-stage hierarchical testing
#' and \code{\link{OTC}} for finding the optimal testing configuration for a
#' number of standard group testing algorithms.
#'
#' \url{http://chrisbilder.com/grouptesting}
#'
#' @family OTC functions
#'
#' @examples
#' # Find the OTC for informative three-stage hierarchical 
#' #   testing over a range of group sizes.
#' # A vector of individual probabilities is generated using
#' #   the expected value of order statistics from a beta 
#' #   distribution with p = 0.05 and a heterogeneity level 
#' #   of alpha = 0.5. Depending on the specified probability, 
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
#' set.seed(8318)
#' Inf.D3(p=0.05, Se=0.99, Sp=0.99, group.sz=3:30, 
#' obj.fn=c("ET", "MAR"), alpha=0.5)}
#' 
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' set.seed(8318)
#' Inf.D3(p=0.05, Se=0.99, Sp=0.99, group.sz=10:15, 
#' obj.fn=c("ET", "MAR"), alpha=0.5)
#'
#' # Find the OTC out of all possible testing configurations
#' #   for a specified group size and vector of individual 
#' #   probabilities.
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' set.seed(1216)
#' p.vec <- p.vec.func(p=0.10, alpha=2, grp.sz=12)
#' Inf.D3(p=p.vec, Se=0.99, Sp=0.99, group.sz=12,
#' obj.fn=c("ET", "MAR", "GR"), weights=matrix(data=c(1,1), 
#' nrow=1, ncol=2, byrow=TRUE), alpha=NA)

# Brianna Hitt - 05-01-17
# Updated: Brianna Hitt - 06-20-18
Inf.D3 <- function(p, Se, Sp, group.sz, obj.fn, weights=NULL, alpha=2){
  
  start.time<-proc.time()
  
  set.of.I <- group.sz
  
  save.ET <- matrix(data=NA, nrow=length(set.of.I), ncol=2*max(set.of.I)+10)
  save.MAR <- matrix(data=NA, nrow=length(set.of.I), ncol=2*max(set.of.I)+10)
  save.GR1 <- matrix(data=NA, nrow=length(set.of.I), ncol=2*max(set.of.I)+10)
  save.GR2 <- matrix(data=NA, nrow=length(set.of.I), ncol=2*max(set.of.I)+10)
  save.GR3 <- matrix(data=NA, nrow=length(set.of.I), ncol=2*max(set.of.I)+10)
  save.GR4 <- matrix(data=NA, nrow=length(set.of.I), ncol=2*max(set.of.I)+10)
  save.GR5 <- matrix(data=NA, nrow=length(set.of.I), ncol=2*max(set.of.I)+10)
  save.GR6 <- matrix(data=NA, nrow=length(set.of.I), ncol=2*max(set.of.I)+10)
  
  count <- 1
  
  for(I in set.of.I){
    # build a vector of probabilities for a heterogeneous population
    if(length(p)==1){
      p.vec <- p.vec.func(p=p, alpha=alpha, grp.sz=I)
    } else if(length(p)>1){
      p.vec <- sort(p)
      alpha <- NA
    }
    
    # generate a matrix of all possible configurations/sets of pool sizes
    # the parts() function requires loading the partitions library
    # do not include the first column, which would test the initial group twice
    possible.groups <- parts(n=I)[,-1]
    
    save.it <- matrix(data=NA, nrow=ncol(possible.groups), ncol=2*max(set.of.I)+17)
    
    counter <- 1
    for(c in 1:ncol(possible.groups)){
      # extract the configuration, ordering, and group sizes for each column
      config <- c(possible.groups[,c], 1:I)
      order.for.p <- config[(1+I):(2*I)]
      gp.sizes <- config[1:I]
      
      # call hierarchical.desc2() for the configuration
      save.info <- hierarchical.desc2(p=p.vec[order.for.p], se=Se, sp=Sp, I2=gp.sizes[gp.sizes!=0], order.p=FALSE)
      
      # extract accuracy measures for each individual
      ET <- save.info$ET
      PSe.vec <- save.info$individual.testerror$pse.vec
      PSp.vec <- save.info$individual.testerror$psp.vec
      if("MAR" %in% obj.fn){
        MAR <- MAR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec)
      } else{MAR <- NA}
      
      # calculate overall accuracy measures
      group.testerror <- save.info$group.testerror
      names(group.testerror) <- NULL
      PSe <- group.testerror[1]
      PSp <- group.testerror[2]
      PPPV <- group.testerror[3]
      PNPV <- group.testerror[4]
      
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
      
      save.it[counter,] <- c(p.vec, rep(NA, max(0, max(set.of.I)-length(p.vec))), alpha, Se, Sp, I, ET, ET/I, MAR, GR1/I, GR2/I, GR3/I, GR4/I, GR5/I, GR6/I, PSe, PSp, PPPV, PNPV, gp.sizes, rep(0, max(0, max(set.of.I)-length(gp.sizes))))
      counter <- counter + 1
      
    }
    
    # find the best configuration for each initial group size I, out of all possible configurations
    save.ET[count,] <- save.it[save.it[,(max(set.of.I)+6)]==min(save.it[,(max(set.of.I)+6)]), c(1:(max(set.of.I)+5),(max(set.of.I)+6),(max(set.of.I)+14):ncol(save.it))]
    if(class(try(save.MAR[count,] <- save.it[save.it[,(max(set.of.I)+7)]==min(save.it[,(max(set.of.I)+7)]), c(1:(max(set.of.I)+5),(max(set.of.I)+7),(max(set.of.I)+14):ncol(save.it))],silent=T))!="try-error"){
      save.MAR[count,] <- save.it[save.it[,(max(set.of.I)+7)]==min(save.it[,(max(set.of.I)+7)]), c(1:(max(set.of.I)+5),(max(set.of.I)+7),(max(set.of.I)+14):ncol(save.it))]
    }
    if(class(try(save.GR1[count,] <- save.it[save.it[,(max(set.of.I)+8)]==min(save.it[,(max(set.of.I)+8)]), c(1:(max(set.of.I)+5),(max(set.of.I)+8),(max(set.of.I)+14):ncol(save.it))],silent=T))!="try-error"){
      save.GR1[count,] <- save.it[save.it[,(max(set.of.I)+8)]==min(save.it[,(max(set.of.I)+8)]), c(1:(max(set.of.I)+5),(max(set.of.I)+8),(max(set.of.I)+14):ncol(save.it))]
    }
    if(class(try(save.GR2[count,] <- save.it[save.it[,(max(set.of.I)+9)]==min(save.it[,(max(set.of.I)+9)]), c(1:(max(set.of.I)+5),(max(set.of.I)+9),(max(set.of.I)+14):ncol(save.it))],silent=T))!="try-error"){
      save.GR2[count,] <- save.it[save.it[,(max(set.of.I)+9)]==min(save.it[,(max(set.of.I)+9)]), c(1:(max(set.of.I)+5),(max(set.of.I)+9),(max(set.of.I)+14):ncol(save.it))]
    }
    if(class(try(save.GR3[count,] <- save.it[save.it[,(max(set.of.I)+10)]==min(save.it[,(max(set.of.I)+10)]), c(1:(max(set.of.I)+5),(max(set.of.I)+10),(max(set.of.I)+14):ncol(save.it))],silent=T))!="try-error"){
      save.GR3[count,] <- save.it[save.it[,(max(set.of.I)+10)]==min(save.it[,(max(set.of.I)+10)]), c(1:(max(set.of.I)+5),(max(set.of.I)+10),(max(set.of.I)+14):ncol(save.it))]
    }
    if(class(try(save.GR4[count,] <- save.it[save.it[,(max(set.of.I)+11)]==min(save.it[,(max(set.of.I)+11)]), c(1:(max(set.of.I)+5),(max(set.of.I)+11),(max(set.of.I)+14):ncol(save.it))],silent=T))!="try-error"){
      save.GR4[count,] <- save.it[save.it[,(max(set.of.I)+11)]==min(save.it[,(max(set.of.I)+11)]), c(1:(max(set.of.I)+5),(max(set.of.I)+11),(max(set.of.I)+14):ncol(save.it))]
    }
    if(class(try( save.GR5[count,] <- save.it[save.it[,(max(set.of.I)+12)]==min(save.it[,(max(set.of.I)+12)]), c(1:(max(set.of.I)+5),(max(set.of.I)+12),(max(set.of.I)+14):ncol(save.it))],silent=T))!="try-error"){
      save.GR5[count,] <- save.it[save.it[,(max(set.of.I)+12)]==min(save.it[,(max(set.of.I)+12)]), c(1:(max(set.of.I)+5),(max(set.of.I)+12),(max(set.of.I)+14):ncol(save.it))]
    }
    if(class(try(save.GR6[count,] <- save.it[save.it[,(max(set.of.I)+13)]==min(save.it[,(max(set.of.I)+13)]), c(1:(max(set.of.I)+5),(max(set.of.I)+13),(max(set.of.I)+14):ncol(save.it))],silent=T))!="try-error"){
      save.GR6[count,] <- save.it[save.it[,(max(set.of.I)+13)]==min(save.it[,(max(set.of.I)+13)]), c(1:(max(set.of.I)+5),(max(set.of.I)+13),(max(set.of.I)+14):ncol(save.it))]    
    }
    
    cat("Initial Group Size =", I, "\n")
    count <- count + 1
  }
  
  # find the optimal testing configuration, over all initial group sizes considered
  result.ET <- save.ET[save.ET[,(max(set.of.I)+6)]==min(save.ET[,(max(set.of.I)+6)]),]
  result.MAR <- save.MAR[save.MAR[,(max(set.of.I)+6)]==min(save.MAR[,(max(set.of.I)+6)]),]
  result.GR1 <- save.GR1[save.GR1[,(max(set.of.I)+6)]==min(save.GR1[,(max(set.of.I)+6)]),]
  result.GR2 <- save.GR2[save.GR2[,(max(set.of.I)+6)]==min(save.GR2[,(max(set.of.I)+6)]),]
  result.GR3 <- save.GR3[save.GR3[,(max(set.of.I)+6)]==min(save.GR3[,(max(set.of.I)+6)]),]
  result.GR4 <- save.GR4[save.GR4[,(max(set.of.I)+6)]==min(save.GR4[,(max(set.of.I)+6)]),]
  result.GR5 <- save.GR5[save.GR5[,(max(set.of.I)+6)]==min(save.GR5[,(max(set.of.I)+6)]),]
  result.GR6 <- save.GR6[save.GR6[,(max(set.of.I)+6)]==min(save.GR6[,(max(set.of.I)+6)]),]
  
  p.vec.ET <- (result.ET[1:max(set.of.I)])[!is.na(result.ET[1:max(set.of.I)])]
  if("MAR" %in% obj.fn){
    p.vec.MAR <- (result.MAR[1:max(set.of.I)])[!is.na(result.MAR[1:max(set.of.I)])]
  } else{p.vec.MAR <- NA}
  if(is.null(dim(weights))){
    p.vec.GR1 <- NA
    p.vec.GR2 <- NA
    p.vec.GR3 <- NA
    p.vec.GR4 <- NA
    p.vec.GR5 <- NA
    p.vec.GR6 <- NA
  } else{
    p.vec.GR1 <- (result.GR1[1:max(set.of.I)])[!is.na(result.GR1[1:max(set.of.I)])]
    if(dim(weights)[1]>=2){
      p.vec.GR2 <- (result.GR2[1:max(set.of.I)])[!is.na(result.GR2[1:max(set.of.I)])]
    } else{p.vec.GR2 <- NA}
    if(dim(weights)[1]>=3){
      p.vec.GR3 <- (result.GR3[1:max(set.of.I)])[!is.na(result.GR3[1:max(set.of.I)])]
    } else{p.vec.GR3 <- NA}
    if(dim(weights)[1]>=4){
      p.vec.GR4 <- (result.GR4[1:max(set.of.I)])[!is.na(result.GR4[1:max(set.of.I)])]
    } else{p.vec.GR4 <- NA}
    if(dim(weights)[1]>=5){
      p.vec.GR5 <- (result.GR5[1:max(set.of.I)])[!is.na(result.GR5[1:max(set.of.I)])]
    } else{p.vec.GR5 <- NA}
    if(dim(weights)[1]>=6){
      p.vec.GR6 <- (result.GR6[1:max(set.of.I)])[!is.na(result.GR6[1:max(set.of.I)])]
    } else{p.vec.GR6 <- NA}
  }
  
  # create a list of results for each objective function
  opt.ET <- list("OTC"=list("Stage1"=result.ET[(max(set.of.I)+4)], "Stage2"=(result.ET[(max(set.of.I)+11):length(result.ET)])[result.ET[(max(set.of.I)+11):length(result.ET)]!=0]), "p.vec"=p.vec.ET,
                 "ET"=result.ET[(max(set.of.I)+5)], "value"=result.ET[(max(set.of.I)+6)], "PSe"=result.ET[(max(set.of.I)+7)], "PSp"=result.ET[(max(set.of.I)+8)], "PPPV"=result.ET[(max(set.of.I)+9)], "PNPV"=result.ET[(max(set.of.I)+10)])
  opt.MAR <- list("OTC"=list("Stage1"=result.MAR[(max(set.of.I)+4)], "Stage2"=(result.MAR[(max(set.of.I)+11):length(result.MAR)])[result.MAR[(max(set.of.I)+11):length(result.MAR)]!=0]), "p.vec"=p.vec.MAR,
                  "ET"=result.MAR[(max(set.of.I)+5)], "value"=result.MAR[(max(set.of.I)+6)], "PSe"=result.MAR[(max(set.of.I)+7)], "PSp"=result.MAR[(max(set.of.I)+8)], "PPPV"=result.MAR[(max(set.of.I)+9)], "PNPV"=result.MAR[(max(set.of.I)+10)])  
  opt.GR1 <- list("OTC"=list("Stage1"=result.GR1[(max(set.of.I)+4)], "Stage2"=(result.GR1[(max(set.of.I)+11):length(result.GR1)])[result.GR1[(max(set.of.I)+11):length(result.GR1)]!=0]), "p.vec"=p.vec.GR1,
                  "ET"=result.GR1[(max(set.of.I)+5)], "value"=result.GR1[(max(set.of.I)+6)], "PSe"=result.GR1[(max(set.of.I)+7)], "PSp"=result.GR1[(max(set.of.I)+8)], "PPPV"=result.GR1[(max(set.of.I)+9)], "PNPV"=result.GR1[(max(set.of.I)+10)])
  opt.GR2 <- list("OTC"=list("Stage1"=result.GR2[(max(set.of.I)+4)], "Stage2"=(result.GR2[(max(set.of.I)+11):length(result.GR2)])[result.GR2[(max(set.of.I)+11):length(result.GR2)]!=0]), "p.vec"=p.vec.GR2,
                  "ET"=result.GR2[(max(set.of.I)+5)], "value"=result.GR2[(max(set.of.I)+6)], "PSe"=result.GR2[(max(set.of.I)+7)], "PSp"=result.GR2[(max(set.of.I)+8)], "PPPV"=result.GR2[(max(set.of.I)+9)], "PNPV"=result.GR2[(max(set.of.I)+10)])  
  opt.GR3 <- list("OTC"=list("Stage1"=result.GR3[(max(set.of.I)+4)], "Stage2"=(result.GR3[(max(set.of.I)+11):length(result.GR3)])[result.GR3[(max(set.of.I)+11):length(result.GR3)]!=0]), "p.vec"=p.vec.GR3,
                  "ET"=result.GR3[(max(set.of.I)+5)], "value"=result.GR3[(max(set.of.I)+6)], "PSe"=result.GR3[(max(set.of.I)+7)], "PSp"=result.GR3[(max(set.of.I)+8)], "PPPV"=result.GR3[(max(set.of.I)+9)], "PNPV"=result.GR3[(max(set.of.I)+10)])
  opt.GR4 <- list("OTC"=list("Stage1"=result.GR4[(max(set.of.I)+4)], "Stage2"=(result.GR4[(max(set.of.I)+11):length(result.GR4)])[result.GR4[(max(set.of.I)+11):length(result.GR4)]!=0]), "p.vec"=p.vec.GR4,
                  "ET"=result.GR4[(max(set.of.I)+5)], "value"=result.GR4[(max(set.of.I)+6)], "PSe"=result.GR4[(max(set.of.I)+7)], "PSp"=result.GR4[(max(set.of.I)+8)], "PPPV"=result.GR4[(max(set.of.I)+9)], "PNPV"=result.GR4[(max(set.of.I)+10)])
  opt.GR5 <- list("OTC"=list("Stage1"=result.GR5[(max(set.of.I)+4)], "Stage2"=(result.GR5[(max(set.of.I)+11):length(result.GR5)])[result.GR5[(max(set.of.I)+11):length(result.GR5)]!=0]), "p.vec"=p.vec.GR5,
                  "ET"=result.GR5[(max(set.of.I)+5)], "value"=result.GR5[(max(set.of.I)+6)], "PSe"=result.GR5[(max(set.of.I)+7)], "PSp"=result.GR5[(max(set.of.I)+8)], "PPPV"=result.GR5[(max(set.of.I)+9)], "PNPV"=result.GR5[(max(set.of.I)+10)])
  opt.GR6 <- list("OTC"=list("Stage1"=result.GR6[(max(set.of.I)+4)], "Stage2"=(result.GR6[(max(set.of.I)+11):length(result.GR6)])[result.GR6[(max(set.of.I)+11):length(result.GR6)]!=0]), "p.vec"=p.vec.GR6,
                  "ET"=result.GR6[(max(set.of.I)+5)], "value"=result.GR6[(max(set.of.I)+6)], "PSe"=result.GR6[(max(set.of.I)+7)], "PSp"=result.GR6[(max(set.of.I)+8)], "PPPV"=result.GR6[(max(set.of.I)+9)], "PNPV"=result.GR6[(max(set.of.I)+10)])
  
  # create a list of results, including all objective functions
  opt.all <- list("opt.ET"=opt.ET, "opt.MAR"=opt.MAR, "opt.GR1"=opt.GR1, "opt.GR2"=opt.GR2, 
                  "opt.GR3"=opt.GR3, "opt.GR4"=opt.GR4, "opt.GR5"=opt.GR5, "opt.GR6"=opt.GR6)
  # remove any objective functions not requested by the user
  opt.req <- Filter(function(x) !is.na(x$ET), opt.all)
  time.it(start.time)
  c("prob"=list(p), "alpha"=alpha, "Se"=Se, "Sp"=Sp, opt.req)
  
}

###################################################################



# Start  NI.D3() function
###################################################################

#' @title Find the optimal testing configuration for non-informative 
#' three-stage hierarchical testing
#'
#' @description Find the optimal testing configuration (OTC) for
#' non-informative three-stage hierarchical testing and calculate
#' the associated operating characteristics.
#'
#' @param p the probability of disease, which can be specified as an overall
#' probability of disease or a homogeneous vector of individual probabilities.
#' @param group.sz a single group size over which to find the OTC
#' out of all possible testing configurations, or a range of
#' group sizes over which to find the OTC.
#' @inheritParams OTC
#'
#' @details This function finds the OTC and computes the associated operating 
#' characteristics for non-informative three-stage hierarchical testing, 
#' by considering all possible testing configurations. 
#' Operating characteristics calculated are expected number of tests, pooling 
#' sensitivity, pooling specificity, pooling positive predictive value, and 
#' pooling negative predictive value for the algorithm. See Hitt et al. (2018) at
#' \url{http://chrisbilder.com/grouptesting} or Kim et al. (2007) 
#' for additional details on the implementation of non-informative three-stage 
#' hierarchical testing.
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
#' \item{Stage1}{pool size for the first stage of testing, i.e. the initial
#' group size.}
#' \item{Stage2}{pool sizes for the second stage of testing.}}}
#' \item{p.vec}{the vector of individual probabilities}
#' \item{ET}{the expected testing expenditure for the OTC.}
#' \item{value}{the value of the objective function per individual.}
#' \item{PSe}{the overall pooling sensitivity for the algorithm. Further 
#' details are given under 'Details'.}
#' \item{PSp}{the overall pooling specificity for the algorithm. Further 
#' details are given under 'Details'.}
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
#' \insertRef{Kim2007}{binGroup}
#' 
#' \insertRef{Malinovsky2016}{binGroup}
#'
#' @seealso
#' \code{\link{Inf.D3}} for informative three-stage hierarchical testing, 
#' \code{\link{NI.Dorf}} for non-informative two-stage hierarchical (Dorfman) 
#' testing, \code{\link{Inf.Dorf}} for informative two-stage hierarchical 
#' testing, and \code{\link{OTC}} for finding the optimal testing configuration 
#' for a number of standard group testing algorithms.
#'
#' \url{http://chrisbilder.com/grouptesting}
#'
#' @family OTC functions
#'
#' @examples
#' # Find the OTC for non-informative three-stage 
#' #   hierarchical testing over a range of group sizes.
#' # This example takes approximately 20 seconds to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' \dontrun{
#' NI.D3(p=0.02, Se=0.90, Sp=0.90, group.sz=3:30, 
#' obj.fn=c("ET", "MAR"))}
#' 
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' NI.D3(p=0.02, Se=0.90, Sp=0.90, group.sz=3:12, 
#' obj.fn=c("ET", "MAR"))
#'
#' # Find the OTC out of all possible configurations for 
#' #   a specified group size for non-informative
#' #   three-stage hierarchical testing.
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' NI.D3(p=rep(0.005, 15), Se=0.99, Sp=0.99, group.sz=15,
#' obj.fn=c("ET", "MAR", "GR"), weights=matrix(data=c(1,1,10,10), 
#' nrow=2, ncol=2, byrow=TRUE))

# Brianna Hitt - 05-01-17
# Updated: Brianna Hitt - 06-20-18
NI.D3 <- function(p, Se, Sp, group.sz, obj.fn, weights=NULL){
  
  start.time<-proc.time()
  
  set.of.I <- group.sz
  
  save.ET <- matrix(data=NA, nrow=length(set.of.I), ncol=max(set.of.I)+10)
  save.MAR <- matrix(data=NA, nrow=length(set.of.I), ncol=max(set.of.I)+10)
  save.GR1 <- matrix(data=NA, nrow=length(set.of.I), ncol=max(set.of.I)+10)
  save.GR2 <- matrix(data=NA, nrow=length(set.of.I), ncol=max(set.of.I)+10)
  save.GR3 <- matrix(data=NA, nrow=length(set.of.I), ncol=max(set.of.I)+10)
  save.GR4 <- matrix(data=NA, nrow=length(set.of.I), ncol=max(set.of.I)+10)
  save.GR5 <- matrix(data=NA, nrow=length(set.of.I), ncol=max(set.of.I)+10)
  save.GR6 <- matrix(data=NA, nrow=length(set.of.I), ncol=max(set.of.I)+10)
  
  count <- 1
  
  for(I in set.of.I){
    
    p.vec <- rep(x=p[1], times=I)
    
    # generate a matrix of all possible configurations/sets of pool sizes
    # the parts() function requires loading the partitions library
    # do not include the first column, which would test the initial group twice
    possible.groups <- parts(n=I)[,-1]
    
    save.it <- matrix(data=NA, nrow=ncol(possible.groups), ncol=max(set.of.I)+17)
    
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
      
      save.it[counter,] <- c(p[1], Se, Sp, I, ET, ET/I, MAR, GR1/I, GR2/I, GR3/I, GR4/I, GR5/I, GR6/I, PSe, PSp, PPPV, PNPV, gp.sizes, rep(0, max(0, max(set.of.I)-length(gp.sizes))))
      counter <- counter + 1
      
    }
    
    # find the best configuration for each initial group size I, out of all possible configurations
    save.ET[count,] <- save.it[save.it[,6]==min(save.it[,6]), c(1:5,6,14:ncol(save.it))]
    if(class(try(save.MAR[count,] <- save.it[save.it[,7]==min(save.it[,7]), c(1:5,7,14:ncol(save.it))],silent=T))!="try-error"){
      save.MAR[count,] <- save.it[save.it[,7]==min(save.it[,7]), c(1:5,7,14:ncol(save.it))]
    }
    if(class(try(save.GR1[count,] <- save.it[save.it[,8]==min(save.it[,8]), c(1:5,8,14:ncol(save.it))],silent=T))!="try-error"){
      save.GR1[count,] <- save.it[save.it[,8]==min(save.it[,8]), c(1:5,8,14:ncol(save.it))]
    }
    if(class(try(save.GR2[count,] <- save.it[save.it[,9]==min(save.it[,9]), c(1:5,9,14:ncol(save.it))],silent=T))!="try-error"){
      save.GR2[count,] <- save.it[save.it[,9]==min(save.it[,9]), c(1:5,9,14:ncol(save.it))]
    }
    if(class(try(save.GR3[count,] <- save.it[save.it[,10]==min(save.it[,10]), c(1:5,10,14:ncol(save.it))],silent=T))!="try-error"){
      save.GR3[count,] <- save.it[save.it[,10]==min(save.it[,10]), c(1:5,10,14:ncol(save.it))]
    }
    if(class(try(save.GR4[count,] <- save.it[save.it[,11]==min(save.it[,11]), c(1:5,11,14:ncol(save.it))],silent=T))!="try-error"){
      save.GR4[count,] <- save.it[save.it[,11]==min(save.it[,11]), c(1:5,11,14:ncol(save.it))]
    }
    if(class(try( save.GR5[count,] <- save.it[save.it[,12]==min(save.it[,12]), c(1:5,12,14:ncol(save.it))],silent=T))!="try-error"){
      save.GR5[count,] <- save.it[save.it[,12]==min(save.it[,12]), c(1:5,12,14:ncol(save.it))]
    }
    if(class(try(save.GR6[count,] <- save.it[save.it[,13]==min(save.it[,13]), c(1:5,13,14:ncol(save.it))],silent=T))!="try-error"){
      save.GR6[count,] <- save.it[save.it[,13]==min(save.it[,13]), c(1:5,13,14:ncol(save.it))]    
    }
    
    cat("Initial Group Size =", I, "\n")
    count <- count + 1
  }
  
  # find the best initial group size I, out of all possible group sizes
  result.ET <- save.ET[save.ET[,6]==min(save.ET[,6]),]
  result.MAR <- save.MAR[save.MAR[,6]==min(save.MAR[,6]),]
  result.GR1 <- save.GR1[save.GR1[,6]==min(save.GR1[,6]),]
  result.GR2 <- save.GR2[save.GR2[,6]==min(save.GR2[,6]),]
  result.GR3 <- save.GR3[save.GR3[,6]==min(save.GR3[,6]),]
  result.GR4 <- save.GR4[save.GR4[,6]==min(save.GR4[,6]),]
  result.GR5 <- save.GR5[save.GR5[,6]==min(save.GR5[,6]),]
  result.GR6 <- save.GR6[save.GR6[,6]==min(save.GR6[,6]),]
  
  p.vec.ET <- rep(x=result.ET[1], times=result.ET[4])
  if("MAR" %in% obj.fn){
    p.vec.MAR <- rep(x=result.MAR[1], times=result.MAR[4])
  } else{p.vec.MAR <- NA}
  if(is.null(dim(weights))){
    p.vec.GR1 <- NA
    p.vec.GR2 <- NA
    p.vec.GR3 <- NA
    p.vec.GR4 <- NA
    p.vec.GR5 <- NA
    p.vec.GR6 <- NA
  } else{
    p.vec.GR1 <- rep(x=result.GR1[1], times=result.GR1[4])
    if(dim(weights)[1]>=2){
      p.vec.GR2 <- rep(x=result.GR2[1], times=result.GR2[4])
    } else{p.vec.GR2 <- NA}
    if(dim(weights)[1]>=3){
      p.vec.GR3 <- rep(x=result.GR3[1], times=result.GR3[4])
    } else{p.vec.GR3 <- NA}
    if(dim(weights)[1]>=4){
      p.vec.GR4 <- rep(x=result.GR4[1], times=result.GR4[4])
    } else{p.vec.GR4 <- NA}
    if(dim(weights)[1]>=5){
      p.vec.GR5 <- rep(x=result.GR5[1], times=result.GR5[4])
    } else{p.vec.GR5 <- NA}
    if(dim(weights)[1]>=6){
      p.vec.GR6 <- rep(x=result.GR6[1], times=result.GR6[4])
    } else{p.vec.GR6 <- NA}
  }
  
  # create a list of results for each objective function
  opt.ET <- list("OTC"=list("Stage1"=result.ET[4], "Stage2"=(result.ET[11:length(result.ET)])[result.ET[11:length(result.ET)]!=0]), "p.vec"=p.vec.ET,
                 "ET"=result.ET[5], "value"=result.ET[6], "PSe"=result.ET[7], "PSp"=result.ET[8], "PPPV"=result.ET[9], "PNPV"=result.ET[10])
  opt.MAR <- list("OTC"=list("Stage1"=result.MAR[4], "Stage2"=(result.MAR[11:length(result.MAR)])[result.MAR[11:length(result.MAR)]!=0]), "p.vec"=p.vec.MAR,
                  "ET"=result.MAR[5], "value"=result.MAR[6], "PSe"=result.MAR[7], "PSp"=result.MAR[8], "PPPV"=result.MAR[9], "PNPV"=result.MAR[10])  
  opt.GR1 <- list("OTC"=list("Stage1"=result.GR1[4], "Stage2"=(result.GR1[11:length(result.GR1)])[result.GR1[11:length(result.GR1)]!=0]), "p.vec"=p.vec.GR1,
                  "ET"=result.GR1[5], "value"=result.GR1[6], "PSe"=result.GR1[7], "PSp"=result.GR1[8], "PPPV"=result.GR1[9], "PNPV"=result.GR1[10])
  opt.GR2 <- list("OTC"=list("Stage1"=result.GR2[4], "Stage2"=(result.GR2[11:length(result.GR2)])[result.GR2[11:length(result.GR2)]!=0]), "p.vec"=p.vec.GR2,
                  "ET"=result.GR2[5], "value"=result.GR2[6], "PSe"=result.GR2[7], "PSp"=result.GR2[8], "PPPV"=result.GR2[9], "PNPV"=result.GR2[10])  
  opt.GR3 <- list("OTC"=list("Stage1"=result.GR3[4], "Stage2"=(result.GR3[11:length(result.GR3)])[result.GR3[11:length(result.GR3)]!=0]), "p.vec"=p.vec.GR3,
                  "ET"=result.GR3[5], "value"=result.GR3[6], "PSe"=result.GR3[7], "PSp"=result.GR3[8], "PPPV"=result.GR3[9], "PNPV"=result.GR3[10])
  opt.GR4 <- list("OTC"=list("Stage1"=result.GR4[4], "Stage2"=(result.GR4[11:length(result.GR4)])[result.GR4[11:length(result.GR4)]!=0]), "p.vec"=p.vec.GR4,
                  "ET"=result.GR4[5], "value"=result.GR4[6], "PSe"=result.GR4[7], "PSp"=result.GR4[8], "PPPV"=result.GR4[9], "PNPV"=result.GR4[10])
  opt.GR5 <- list("OTC"=list("Stage1"=result.GR5[4], "Stage2"=(result.GR5[11:length(result.GR5)])[result.GR5[11:length(result.GR5)]!=0]), "p.vec"=p.vec.GR5,
                  "ET"=result.GR5[5], "value"=result.GR5[6], "PSe"=result.GR5[7], "PSp"=result.GR5[8], "PPPV"=result.GR5[9], "PNPV"=result.GR5[10])
  opt.GR6 <- list("OTC"=list("Stage1"=result.GR6[4], "Stage2"=(result.GR6[11:length(result.GR6)])[result.GR6[11:length(result.GR6)]!=0]), "p.vec"=p.vec.GR6,
                  "ET"=result.GR6[5], "value"=result.GR6[6], "PSe"=result.GR6[7], "PSp"=result.GR6[8], "PPPV"=result.GR6[9], "PNPV"=result.GR6[10])
  
  # create a list of results, including all objective functions
  opt.all <- list("opt.ET"=opt.ET, "opt.MAR"=opt.MAR, "opt.GR1"=opt.GR1, "opt.GR2"=opt.GR2, 
                  "opt.GR3"=opt.GR3, "opt.GR4"=opt.GR4, "opt.GR5"=opt.GR5, "opt.GR6"=opt.GR6)
  # remove any objective functions not requested by the user
  opt.req <- Filter(function(x) !is.na(x$ET), opt.all)
  time.it(start.time)
  c("prob"=list(p), "Se"=Se, "Sp"=Sp, opt.req)
  
}

###################################################################


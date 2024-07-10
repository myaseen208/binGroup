#############################################################################
# This file contains functions to calculate the expected testing expenditure 
# of array testing, and the individual specific measures of accuracy Pooling 
# Sensitivity, Pooling Specificity, Pooling Positive Predicitive Value, and 
# Pooling Negative Predicitive Value.
# Last modified date: 8-15-2011
# Author: Chris McMahan
#############################################################################

#' @title Operating characteristics for array testing without master pooling
#' 
#' @description Calculate the expected number of tests and accuracy measures
#' for each individual using array testing without master pooling
#' 
#' @param p matrix of probabilities corresponding to each individual's risk 
#' of disease.
#' @inheritParams hierarchical.desc2
#' 
#' @details This function calculates the operating characteristics for
#' array testing without master pooling. Operating characteristics calculated
#' are expected number of tests, pooling sensitivity, pooling specificity, 
#' pooling positive predictive value, and pooling negative predictive value 
#' for each individual. 
#'
#' @return A list containing:
#' \item{T}{the expected number of tests for the array.}
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
#' @author This function was originally written by Christopher S. McMahan for 
#' McMahan et al. (2012). The function was obtained from 
#' \url{http://chrisbilder.com/grouptesting}.
#' 
#' @references 
#' \insertRef{McMahan2012b}{binGroup}
#' 
#' @seealso
#' \code{\link{MasterPool.Array.Measures}} for calculating operating 
#' characteristics under non-informative array testing with master pooling, 
#' \code{\link{hierarchical.desc2}} for three-stage hierarchical and 
#' non-informative two-stage hierarchical testing, and 
#' \code{\link{inf.dorf.measures}} for informative two-stage hierarchical 
#' testing. See \code{\link{p.vec.func}} for generating a vector of 
#' individual risk probabilities for informative group testing and 
#' \code{\link{Informative.array.prob}} for arranging individual risk
#' probabilities in a matrix for informative array testing.  
#' 
#' \url{http://chrisbilder.com/grouptesting}
#' 
#' @family Operating characteristic functions
#' 
#' @examples 
#' # Calculate the operating characteristics for 
#' #   non-informative array testing without master
#' #   pooling, with a 5x5 array and an overall disease 
#' #   risk of p = 0.02.
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' p.mat <- matrix(data=0.02, ncol=5, nrow=5)
#' Array.Measures(p=p.mat, se=0.95, sp=0.95)
#' 
#' # Calculate the operating characteristics for 
#' #   informative array testing without master
#' #   pooling, with a 3x3 array and an overall disease
#' #   risk of p = 0.03 and alpha = 2.
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' set.seed(8791)
#' p.vec <- p.vec.func(p=0.03, alpha=2, grp.sz=9)
#' p.mat <- Informative.array.prob(prob.vec=p.vec, nr=3, 
#' nc=3, method="gd")
#' Array.Measures(p=p.mat, se=0.99, sp=0.99)

Array.Measures<-function(p,se,sp){
  d<-dim(p)
  J<-d[1]
  K<-d[2]

  #######################################################
  # Finds probability that each individual tests positive
  e.ind<-matrix(-100,nrow=J,ncol=K)
  r<-matrix(0:J,ncol=1,nrow=(J+1))
  RA<-apply(r,1,rc.pos,p=p,id=1)
  c<-matrix(0:K,ncol=1,nrow=(K+1))
  CB<-apply(c,1,rc.pos,p=p,id=2)

  for(j in 1:J){
    for(k in 1:K){
      
      p.temp1<-p
      p.temp1[,k]<-rep(0,J)
      RAJ<-apply(r,1,rc.pos,p=p.temp1,id=1)
      p.c<-prod((1-p[,k]))
      a1<-sum( (1-sp)*(1-se)^r*sp^(J-r)*p.c*RAJ + se*(1-se)^r*sp^(J-r)*((RA-p.c*RAJ)) )
      
      p.temp1<-p
      p.temp1[j,]<-rep(0,K)
      CBK<-apply(c,1,rc.pos,p=p.temp1,id=2)
      p.r<-prod((1-p[j,]))
      a2<-sum( (1-sp)*(1-se)^c*sp^(K-c)*p.r*CBK + se*(1-se)^c*sp^(K-c)*((CB-p.r*CBK)) )
      
      a3<- se^2+(1-se-sp)^2*prod(1-p[,k])*prod(1-p[j,])/(1-p[j,k])+(se*(1-sp)-se^2)*(prod(1-p[j,])+prod(1-p[,k]))
      
      e.ind[j,k]<-(a1+a2+a3)
    }
  }
  T<-sum(e.ind)+J+K

  ###############################################
  # Finds Pooling Sensitivity for each individual
  pse.ind<-matrix(-100,nrow=J,ncol=K)
  c0<-1-(se+(1-se-sp)*apply((1-p),2,prod))
  r0<-1-(se+(1-se-sp)*apply((1-p),1,prod))
  
  for(j in 1:J){
    for(k in 1:K){
      pse.ind[j,k]<-se^3+se^2*(1-se)*(prod(c0[-k])+prod(r0[-j]))
    }
  }

  ###############################################
  # Finds Pooling Specificity for each individual
  psp.ind<-matrix(-100,nrow=J,ncol=K)
  r<-matrix(0:J,ncol=1,nrow=(J+1))
  c<-matrix(0:K,ncol=1,nrow=(K+1))
  
  for(j in 1:J){
    for(k in 1:K){
      p.temp1<-p
      p.temp1[,k]<-rep(0,J)
      p.temp2<-p
      p.temp2[j,k]<-0
      
      RAJ<-apply(r,1,rc.pos,p=p.temp1,id=1)
      RAj<-apply(r,1,rc.pos,p=p.temp2,id=1)
      p.c<-prod((1-p.temp2[,k]))
      a1<-sum((1-sp)*(1-se)^r*sp^(J-r)*p.c*RAJ + se*(1-se)^r*sp^(J-r)*((RAj-p.c*RAJ)))
      
      p.temp1<-p
      p.temp1[j,]<-rep(0,K)
      p.temp2<-p
      p.temp2[j,k]<-0
      
      CBK<-apply(r,1,rc.pos,p=p.temp1,id=2)
      CBk<-apply(r,1,rc.pos,p=p.temp2,id=2)
      p.r<-prod((1-p.temp2[j,]))
      a2<-sum((1-sp)*(1-se)^c*sp^(K-c)*p.r*CBK + se*(1-se)^c*sp^(K-c)*((CBk-p.r*CBK)))
      
      a3<-(se+(1-se-sp)*prod((1-p[,k]))/(1-p[j,k]))*(se+(1-se-sp)*prod((1-p[j,]))/(1-p[j,k])) 
      psp.ind[j,k]<-1-(1-sp)*(a1+a2+a3)
    }
  }
  
  #########################################################################
  # Finds Pooling Positive (Negative) Predicitive Value for each individual
  ppv.ind<-p*pse.ind/(p*pse.ind+(1-p)*(1-psp.ind))
  npv.ind<-(1-p)*psp.ind/((1-p)*psp.ind + p*(1-pse.ind))
  
  return(list("T"=T,"PSe"=pse.ind,"PSp"=psp.ind, "PPV"=ppv.ind, "NPV"=npv.ind))
}

##################################################################
##################################################################
# Support Functions needed for Array.Measures()               ####
##################################################################
##################################################################
rc.pos<-function(n,p,id){
d<-dim(p)
M<-d[id]
p.pos<-1-apply((1-p),id,prod)
p.frac<-p.pos/(1-p.pos)
  if(n==0){
  res<-prod(1-p.pos) 
  }
  if(n>0){
    if(n==1){
    res<-sum(p.frac)*prod(1-p.pos)
    }
    if(n>1){
    temp<-cumsum(p.frac[M:n])
    temp<-temp[(M-n+1):1]
      if(n>2){
        for(i in 1:(n-2)){
        temp<-p.frac[(n-i):(M-i)]*temp
        temp<-cumsum(temp[(M-n+1):1])
        temp<-temp[(M-n+1):1]
        }
      }
    temp<-sum(p.frac[1:(M-n+1)]*temp)
    res<-temp*prod(1-p.pos)
    }
  }
return(res)
}


###############################################################################
###############################################################################
###### Spiral or Gradient Array Construction Function                 #########
###############################################################################
###############################################################################

#' @title Arrange a matrix of probabilities for informative array testing
#' 
#' @description Arrange a vector of individual risk probabilities in a
#' matrix for informative array testing without master pooling. 
#' 
#' @param prob.vec vector of individual risk probabilities, of 
#' length nr*nc.
#' @param nr the number of rows in the array.
#' @param nc the number of columns in the array.
#' @param method character string defining the method to be used
#' for matrix arrangement. Options include spiral ("\kbd{sd}")
#' and gradient ("\kbd{gd}") arrangement. See McMahan et al. (2012)
#' for additional details.
#' 
#' @return A matrix of probabilities arranged according to the 
#' specified method.
#' 
#' @author This function was originally written by Christopher S. McMahan 
#' for McMahan et al. (2012). The function was obtained from 
#' \url{http://chrisbilder.com/grouptesting}.
#' 
#' @references 
#' \insertRef{McMahan2012b}{binGroup}
#' 
#' @seealso 
#' \code{\link{p.vec.func}} for generating a vector of individual
#' risk probabilities from an overall probability of disease and
#' \code{\link{Array.Measures}} for calculating operating characteristics
#' for informative array testing without master pooling.
#' 
#' \url{http://chrisbilder.com/grouptesting}
#' 
#' @family Individual risk probability functions
#' 
#' @examples
#' # Use the gradient arrangement method to create a matrix
#' #   of individual risk probabilities for a 10x10 array.
#' # Depending on the specified probability, alpha level, 
#' #   and overall group size, simulation may be necessary in 
#' #   order to generate the vector of individual probabilities. 
#' #   This is done using the p.vec.func() function and requires 
#' #   the user to set a seed in order to reproduce results.
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' set.seed(1107)
#' p.vec1 <- p.vec.func(p=0.05, alpha=2, grp.sz=100)
#' Informative.array.prob(prob.vec=p.vec1, nr=10, nc=10, method="gd")
#' 
#' # Use the spiral arrangement method to create a matrix
#' #   of individual risk probabilities for a 5x5 array.
#' # Depending on the specified probability, alpha level, 
#' #   and overall group size, simulation may be necessary in 
#' #   order to generate the vector of individual probabilities. 
#' #   This is done using the p.vec.func() function and requires 
#' #   the user to set a seed in order to reproduce results.
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' set.seed(8791)
#' p.vec2 <- p.vec.func(p=0.02, alpha=0.5, grp.sz=25)
#' Informative.array.prob(prob.vec=p.vec2, nr=5, nc=5, method="sd")

Informative.array.prob<-function(prob.vec,nr,nc, method="sd"){

prob.vec<-sort(prob.vec,decreasing=TRUE)
  if(method=="sd"){
    array.probs<-prob.vec[1:2]
    prob.vec<-prob.vec[-(1:2)]
    array.probs<-cbind(array.probs,sort(prob.vec[1:2],decreasing=FALSE))
    prob.vec<-prob.vec[-(1:2)]
     
    max.iter<-max(nr,nc)

      for(i in 1:max.iter){
        if(nrow(array.probs) < nr){
          array.probs<-rbind(array.probs,prob.vec[1:(ncol(array.probs))])
          prob.vec<-prob.vec[-(1:(ncol(array.probs)))]
        }
        if(ncol(array.probs) < nc){
          array.probs<-cbind(array.probs,sort(prob.vec[1:(nrow(array.probs))],decreasing=FALSE))
          prob.vec<-prob.vec[-(1:(nrow(array.probs)))]
        }
      }  
    }

  if(method=="gd"){
    array.probs<-matrix(prob.vec,ncol=max(nr,nc), nrow=min(nr,nc), byrow=FALSE)
  }
return(array.probs)
}

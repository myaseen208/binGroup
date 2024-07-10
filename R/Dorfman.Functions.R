#######################################################################################
# Used to find the characteristics of an Informative Dorfman decoding process specified  
# by method= OD (Optimal Dorfman), TOD (Thresholded Optimal Dorfman), and PSOD (Pool Specific Optimal Dorfman)
# given p (a vector of all subjects' infection probabilities), 
# se (sensitivity of diagnostic test),
# sp (specificity of diagnostic test), 
# max.pool (maximum allowable pool size),
# thresh.pool (initial pool size used for TOD if threshold is not specified), and 
# threshold (threshold value for TOD, note if a threshold value is not specified one is found algorithmically). 
# Function returns 
# p.star (threshold value used, only applies when using TOD),
# res.e (the expected expenditure of the decoding process),
# res.v (the variance of expenditure of the decoding process), and 
# res.mat a matrix of summary measures which includes each subject's infection probability , pool (pool to which they belong),
# PSe (pooling sensitivity), PSp (pooling specificity), PPV (pooling positive predictive value), and 
# NPV (pooling negative predictive value).
#
#EXAMPLE: opt.info.dorf(prob=rbeta(1000,1,10), se = 1, sp = 1, method ="OD", max.pool=15, thresh.pool=8, threshold=NULL)

#' @title Find the characteristics of an informative two-stage 
#' hierarchical (Dorfman) algorithm
#' 
#' @description Find the characteristics of an informative 
#' two-stage hierarchical (Dorfman) decoding process using Optimal 
#' Dorfman (OD), Thresholded Optimal Dorfman (TOD), or Pool-Specific 
#' Optimal Dorfman (PSOD) algorithms.
#' 
#' @param prob a vector of all subjects' infection probabilities.
#' @param se the sensitivity of the diagnostic test.
#' @param sp the specificity of the diagnostic test.
#' @param method character string defining the specific screening
#' procedure for implementation of Dorfman retesting in a 
#' heterogeneous population. Options include Optimal Dorfman 
#' ("\kbd{OD}"), Thresholded Optimal Dorfman ("\kbd{TOD}"), and
#' Pool-Specific Optimal Dorfman ("\kbd{PSOD}"). Further details 
#' are given under 'Details'.
#' @param max.pool the maximum allowable pool size. Further details 
#' are given under 'Details'.
#' @param thresh.pool the initial pool size used for TOD, if 
#' \kbd{threshold} is not specified. Further details are given 
#' under 'Details'.
#' @param threshold the threshold value for TOD. If a threshold
#' value is not specified, one is found algorithmically. Further 
#' details are given under 'Details'.
#' 
#' @details This function finds the characteristics of an informative
#' two-stage hierarchical (Dorfman) decoding process. Characteristics
#' found include the expected expenditure of the decoding process, 
#' the variance of the expenditure of the decoding process, and the 
#' pooling sensitivity, pooling specificity, pooling positive predictive
#' value, and pooling negative predictive value for each individual.
#' Calculations of these characteristics are done using equations
#' presented in McMahan et al. (2012). 
#' 
#' Optimal Dorfman (OD) is an informative Dorfman algorithm in 
#' which the common pool size \eqn{c=c_{opt}} minimizes 
#' \eqn{E(T^(c))}, the expected number of tests needed to decode 
#' all \eqn{N} individuals when pools of size \eqn{c} are used. 
#' 
#' Thresholded Optimal Dorfman (TOD) is an informative Dorfman 
#' algorithm in which all \eqn{N} individuals are partitioned 
#' into two classes, low-risk and high-risk individuals, based 
#' on whether their risk probability falls below or above a 
#' particular threshold value. The threshold can be specified 
#' using the \kbd{threshold} argument or the TOD algorithm can 
#' identify the optimal threshold value. The low-risk individuals
#' are tested using a optimal common pool size,  and high-risk 
#' individuals are tested individually.
#' 
#' Pool-Specific Optimal Dorfman (PSOD) is an informative Dorfman
#' algorithm in which optimal sizes are determined for each pool. 
#' A total of \eqn{N} individuals are tested in pools that minimize
#' the expected number of tests per individual, on a pool-by-pool 
#' basis. If desired, the user can add the constraint of a maximum 
#' allowable pool size, so that each pool will contain no more 
#' than the maximum allowable number of individuals.
#' 
#' All three informative Dorfman procedures described above require
#' individuals to be ordered from smallest to largest probability 
#' of infection. See McMahan et al. (2012) for additional details
#' on the implementation of informative two-stage hierarchical
#' (Dorfman) testing algorithms.
#' 
#' @return A list containing:
#' \item{tv}{the threshold value used for TOD, if applicable.}
#' \item{e}{the expected expenditure of the decoding process.}
#' \item{v}{the variance of the expenditure of the decoding 
#' process.}
#' \item{summary}{a matrix of summary measures that includes
#' each individual's infection probability, pool (pool to which
#' they belong), pooling sensitivity, pooling specificity, 
#' pooling positive predictive value, and pooling negative
#' predictive value.}
#' 
#' @author This function was originally written by Christopher S. 
#' McMahan for McMahan et al. (2012). The function was obtained 
#' from \url{http://chrisbilder.com/grouptesting}.
#' 
#' @references 
#' \insertRef{Dorfman1943}{binGroup}
#' 
#' \insertRef{McMahan2012a}{binGroup}
#' 
#' @seealso 
#' \url{http://chrisbilder.com/grouptesting}
#' 
#' @family Informative Dorfman functions
#' 
#' @examples 
#' # Find the characteristics of an informative
#' #   Dorfman algorithm, using the OD procedure.
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' opt.info.dorf(prob=rbeta(1000,1,10), se=1, sp=1, 
#' method ="OD", max.pool=15, thresh.pool=8, threshold=NULL)
#' 
#' # Find the characteristics of an informative 
#' #   Dorfman algorithm, using the TOD procedure.
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' set.seed(1002)
#' p.vec <- p.vec.func(p=0.01, alpha=2, grp.sz=20)
#' opt.info.dorf(prob=p.vec, se=0.95, sp=0.95, 
#' method="TOD", max.pool=5, threshold=0.015)

opt.info.dorf<-function(prob, se = 1, sp = 1, method ="OD", max.pool=15, thresh.pool=8, threshold=NULL){

# Saves original ordering
ind.order<-order(prob)

# Orders subjects, required under all Informative measures
prob<-sort(prob)

# Determines number of subjects being screened, and sets up vectors for 
# storing summary measures. Also initializes the threshold p.star
N<-length(prob)
pool.id<-rep(-1000,N)
PSe<-rep(-100,N)
PSp<-rep(-100,N)
PPV<-rep(-100,N)
NPV<-rep(-100,N)
p.star<-threshold


# If method is TOD this finds the threshold value and divides the subjects into 
# the high and low risk classes
if(method=="TOD"){
  if(is.null(p.star)){
  p.star<-thresh.val.dorf(p=prob, psz=thresh.pool, se=se, sp=sp)
  }
  if(p.star < 1){
  N.high<-length(prob[prob > p.star])
  N<-N-N.high
  }
}


# Finds optimal pool size to be used with OD or, finds optimal pool size to be used 
# to decode low risk class in TOD
if(method=="OD" | method=="TOD"){
   if(N==0){psz<-NULL}
   if(N > 0){
     res.psz<-opt.pool.size(p=prob[1:N] ,max.p=max.pool, se=se, sp=sp)
     J<-ceiling(N/res.psz)
     rem<-N-(J-1)*res.psz
     if(rem!=0){
     psz<-c(rep(res.psz,(J-1)),rem)
     }  
     if(rem==0){
     psz<-rep(res.psz,J)
     }
  }
  if(method=="TOD"){
   if(p.star < 1){
   psz<-c(psz,rep(1,N.high))
   J<-length(psz)
   N<-N+N.high
   }
  }
}


# Finds pool sizes to be used with PSOD
if(method=="PSOD"){
  psz<-pool.specific.dorf(p=prob , max.p=max.pool , se=se , sp=sp)
  J<-length(psz)
}

# Finds measures pool by pool
psz<-c(psz,0)
lower<-1
upper<-psz[1]
vec.e<-rep(-1000,J)
vec.v<-rep(-1000,J)
for(i in 1:J){
p.pool<-prob[lower:upper]
pool.id[lower:upper]<-rep(i,length(p.pool))

res<-characteristics.pool(p=p.pool,se=se,sp=sp)
vec.e[i]<-res$e
vec.v[i]<-res$v

res.acc<-accuracy.dorf(p=p.pool,se=se,sp=sp)
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

res.mat<-matrix(c(pool.id, prob, PSe, PSp, PPV, NPV),nrow=N,ncol=6,byrow=FALSE, 
                dimnames=list(as.character(1:N) , c("pool","probability","PSe", "PSp", "PPV", "NPV")))

prob<-prob[ind.order]
return(list("tv"=p.star, "e"=res.e, "v"=res.v, "summary"=res.mat))
}



######################################
### DORFMAN DECODING FUNCTIONS #######
######################################

############################################################################
# Function for the expectation and variation in testing expenditure 
# of a pool of size greater than or equal to one
# here p is a vector of all subjects probability of infection,
# se is sensitivity, and
# sp is specificity.
# res.e is the expected expenditure of the pool and
# res.v is the variance of expenditure of the pool.

#' @title Testing expenditure for informative Dorfman testing
#' 
#' @description Calculate the expectation and variation of the testing
#' expenditure of a pool used with informative Dorfman testing.
#' 
#' @param p a vector of each individual's probability of infection.
#' @inheritParams opt.info.dorf
#' 
#' @details This function calculates the expected value and variance
#' of the testing expenditure of a pool of size greater than or equal 
#' to one used with informative Dorfman testing. Calculations of these 
#' measures are done using the equations presented in McMahan et al. (2012). 
#' 
#' @return a list containing:
#' \item{e}{the expected testing expenditure of the pool.}
#' \item{v}{the variation of the testing expenditure of the pool.}
#'  
#' @author This function was originally written by Christopher S. 
#' McMahan for McMahan et al. (2012). The function was obtained 
#' from \url{http://chrisbilder.com/grouptesting}.
#' 
#' @references 
#' \insertRef{McMahan2012a}{binGroup}
#' 
#' @seealso 
#' \url{http://chrisbilder.com/grouptesting} 
#' 
#' @family Informative Dorfman functions
#' 
#' @examples 
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' set.seed(8135)
#' p.vec <- p.vec.func(p=0.02, alpha=1, grp.sz=10)
#' characteristics.pool(p=p.vec[1:3], se=0.90, sp=0.90)

characteristics.pool<-function(p,se,sp){
n<-length(p)
if(n>1){
prob<-se+(1-se-sp)*prod(1-p)
res.e<-1+n*prob
res.v<-n^2*prob*(1-prob)
}
if(n==1){
res.e<-1
res.v<-0
}
return(list("e"=res.e, "v"=res.v))
}



############################################################################
# Function that returns PSe, PSp, PPV, and NPV for all individuals
# belonging to a pool of size greater than or equal to one
# here p is a vector of all subject probabilities,
# se is sensitivity, and
# sp is specificity.

#' @title Accuracy measures for informative Dorfman testing
#' 
#' @description Calculate the accuracy measures for each individual 
#' in a pool used with informative Dorfman testing.
#' 
#' @param p a vector of each individual's probability of infection.
#' @inheritParams opt.info.dorf
#' 
#' @details This function calculates the pooling sensitivity, pooling
#' specificity, pooling positive predictive value, and pooling negative
#' predictive value for each individual belonging to a pool of size 
#' greater than or equal to one used with informative Dorfman testing. 
#' Calculations of these measures are done using the equations presented 
#' in McMahan et al. (2012). 
#' 
#' @return a list containing:
#' \item{PSe}{a vector containing each individual's pooling sensitivity.}
#' \item{PSp}{a vector containing each individual's pooling specificity.}
#' \item{PPV}{a vector containing each individual's pooling positive 
#' predictive value.}
#' \item{NPV}{a vector containing each individual's pooling negative 
#' predictive value.}
#'  
#' @author This function was originally written by Christopher S. 
#' McMahan for McMahan et al. (2012). The function was obtained 
#' from \url{http://chrisbilder.com/grouptesting}.
#' 
#' @references 
#' \insertRef{McMahan2012a}{binGroup}
#' 
#' @seealso 
#' \url{http://chrisbilder.com/grouptesting}
#' 
#' @family Informative Dorfman functions
#' 
#' @examples 
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' set.seed(8135)
#' p.vec <- p.vec.func(p=0.02, alpha=1, grp.sz=10)
#' accuracy.dorf(p=p.vec[1:3], se=0.90, sp=0.90)

accuracy.dorf<-function(p,se,sp){
cj<-length(p)
se.vec<-rep(-100,cj)
sp.vec<-rep(-100,cj)
ppv.vec<-rep(-100,cj)
npv.vec<-rep(-100,cj)

if(cj==1){
se.vec<-se
sp.vec<-sp
ppv.vec<-p*se/(p*se+(1-p)*(1-sp))
npv.vec<-(1-p)*sp/((1-p)*sp + p*(1-se))
}

if(cj>1){
for(i in 1:cj){
se.vec[i]<-se^2
sp.vec[i]<-1-(1-sp)*(se+(1-se-sp)*prod(1-p[-i]))
ppv.vec[i]<-p[i]*se^2/(p[i]*se^2+(1-p[i])*(1-sp.vec[i]))
npv.vec[i]<-(1-p[i])*sp.vec[i]/((1-p[i])*sp.vec[i] + p[i]*(1-se^2))
}
}

return(list("PSe"=se.vec, "PSp"=sp.vec, "PPV"=ppv.vec, "NPV"=npv.vec))
}



##############################################
# Pooling Algorithm Minimizes Individual 
# Risk on a per pool basis, so it allows for multiple 
# pooling sizes (psz) here p is a vector of all subject probabilities,
# se is sensitivity,
# sp is specificity, and
# max.p is the maximum allowable pool size.

#' @title Find the optimal pool sizes for Pool-Specific Optimal Dorfman 
#' (PSOD) testing
#' 
#' @description Find the set of optimal pool sizes for Pool-Specific 
#' Optimal Dorfman (PSOD) testing.
#' 
#' @param p a vector of each individual's probability of infection.
#' @param max.p the maximum allowable pool size.
#' @inheritParams opt.info.dorf
#' 
#' @details This function finds the set of optimal pool sizes for PSOD
#' testing. PSOD testing uses a greedy algorithm and does not consider
#' all possible sets of pool sizes. See McMahan et al. (2012) for 
#' additional details on the implementation of PSOD testing.
#' 
#' @return The optimal set of pool sizes for PSOD testing.
#'  
#' @author This function was originally written by Christopher S. 
#' McMahan for McMahan et al. (2012). The function was obtained 
#' from \url{http://chrisbilder.com/grouptesting}.
#' 
#' @references 
#' \insertRef{McMahan2012a}{binGroup}
#' 
#' @seealso 
#' \url{http://chrisbilder.com/grouptesting}
#' 
#' @family Informative Dorfman functions
#' 
#' @examples 
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' set.seed(8135)
#' p.vec <- p.vec.func(p=0.02, alpha=1, grp.sz=10)
#' pool.specific.dorf(p=p.vec, max.p=3, se=0.95, sp=0.95)

pool.specific.dorf<-function(p,max.p,se,sp){
p<-sort(p)
N<-length(p)
psz<-rep(-99,N)
k<-0
ind<-1

while(k!=N){
et<-1
max=min(length(p),max.p)
if(length(p)>1){
et<-c(et,rep(99,max))
  for(i in 2:max){
  et[i]<-((i)*(se-(se+sp-1)*prod(1-p[1:(i)]))+1)/(i)
  }}
  m<-1:max
  m<-m[order(et)[1]]
  psz[ind]<-m
  ind<-ind+1
  p<-p[-(1:m)]
  k<-k+m
}
psz<-psz[psz!=-99]
return(psz)
}



####################################################
# Function for finding the optimal pool size (psz),                                               
# here p is a vector of all subject probabilities,
# se is sensitivity,
# sp is specificity, and
# max.p is the maximum allowable pool size.

#' @title Find the optimal pool size for Optimal Dorfman or
#' Thresholded Optimal Dorfman
#' 
#' @description Find the optimal common pool size for Optimal Dorfman
#' (OD) or Thresholded Optimal Dorfman (TOD) testing.
#' 
#' @param p a vector of each individual's probability of infection.
#' @param max.p the maximum allowable pool size.
#' @inheritParams opt.info.dorf
#' 
#' @details This function finds the optimal common pool size for OD or
#' TOD testing. Using OD testing, all individuals are tested using an 
#' optimal common pool size. Using TOD testing, individuals are partitioned 
#' into low-risk and high-risk groups, and all low-risk individuals are 
#' tested using an optimal common pool size. See McMahan et al. (2012) for 
#' additional details on the implementation of OD or TOD testing.
#' 
#' @return The optimal common pool size for OD or TOD testing.
#'  
#' @author This function was originally written by Christopher S. 
#' McMahan for McMahan et al. (2012). The function was obtained 
#' from \url{http://chrisbilder.com/grouptesting}.
#' 
#' @references 
#' \insertRef{McMahan2012a}{binGroup}
#' 
#' @seealso 
#' \url{http://chrisbilder.com/grouptesting}
#' 
#' @family Informative Dorfman functions
#' 
#' @examples 
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' set.seed(8135)
#' p.vec <- p.vec.func(p=0.02, alpha=1, grp.sz=10)
#' opt.pool.size(p=p.vec, max.p=3, se=0.95, sp=0.95)

opt.pool.size<-function(p, max.p, se=1, sp=1){
               
N<-length(p)
M<-min(N, max.p)
p<-sort(p)
e0<-N
psz<-2
L<-ceiling(N/psz)
if(N<(L*psz)){psz.vec<-c(rep(psz,L-1),(N-psz*(L-1)))}
if(N==(L*psz)){psz.vec<-rep(psz,L)}
e1<-0
sub.ind<-c(0,cumsum(psz.vec))
for(m in 1:L){
e1<-e1+(1+psz.vec[m]*(se+(1-se-sp)*prod(1-p[(sub.ind[m]+1):sub.ind[m+1]])))
}

while(e1<e0 & psz<=max.p){
e0<-e1
psz<-psz+1
L<-ceiling(N/psz)
if(N<(L*psz)){psz.vec<-c(rep(psz,L-1),(N-psz*(L-1)))}
if(N==(L*psz)){psz.vec<-rep(psz,L)}
e1<-0
sub.ind<-c(0,cumsum(psz.vec))
for(m in 1:L){
e1<-e1+(1+psz.vec[m]*(se+(1-se-sp)*prod(1-p[(sub.ind[m]+1):sub.ind[m+1]])))
}}

return(psz-1)
}



######################################################################
# Thresholding function for Dorfman Procedure, finds threshold (p.star).                                                                   
# Here p is a vector of all subject probabilities,
# se is sensitivity,
# sp is specificity, and
# psz is the initial pool size.

#' @title Find the optimal threshold value for Thresholded Optimal 
#' Dorfman testing
#' 
#' @description Find the optimal threshold value for Thresholded Optimal
#' Dorfman (TOD) testing.
#' 
#' @param p a vector of each individual's probability of infection.
#' @param psz the initial pool size.
#' @inheritParams opt.info.dorf
#' 
#' @details This function finds the optimal threshold value for TOD 
#' testing for situations where the threshold value is not specified. 
#' See McMahan et al. (2012) for additional details on the implementation 
#' of TOD testing.
#' 
#' @return The optimal threshold value for TOD testing.
#'  
#' @author This function was originally written by Christopher S. 
#' McMahan for McMahan et al. (2012). The function was obtained 
#' from \url{http://chrisbilder.com/grouptesting}.
#' 
#' @references 
#' \insertRef{McMahan2012a}{binGroup}
#' 
#' @seealso 
#' \url{http://chrisbilder.com/grouptesting}
#' 
#' @family Informative Dorfman functions
#' 
#' @examples
#' # This example takes approximately 4 seconds to run. 
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' \dontrun{
#' set.seed(3154)
#' p.vec <- p.vec.func(p=0.10, alpha=0.5, grp.sz=1000)
#' thresh.val.dorf(p=p.vec, psz=10, se=0.95, sp=0.95)}
#' 
#' # This example takes less than 1 second to run.
#' # Estimated running time was calculated using a 
#' #   computer with 16 GB of RAM and one core of an 
#' #   Intel i7-6500U processor.
#' p.vec <- p.vec.func(p=0.15, alpha=2, grp.sz=100)
#' thresh.val.dorf(p=p.vec, psz=10, se=0.95, sp=0.95)

thresh.val.dorf<-function(p, psz, se=1, sp=1){
  N<-length(p)
  p<-sort(p)
  L<-ceiling(N/psz)
  p.star<-1
  upper<-N
  lower<-upper-psz+1
  lower<-max(1,lower)
  Ex.pool<-length(p[lower:upper])*(se+(1-se-sp)*prod(1-p[lower:upper]))+1
    if(Ex.pool > length(p[lower:upper])){
        while(Ex.pool > length(p[lower:upper]) & lower > 1){
        upper<-upper-psz
        lower<-lower-psz
        lower<-max(1,lower)
        Ex.pool<-length(p[lower:upper])*(se+(1-se-sp)*prod(1-p[lower:upper]))+1
        }
    p.star<-(p[upper]+p[upper+1])/2
    }
    if(lower==1){p.star<-0}
return(p.star)

}




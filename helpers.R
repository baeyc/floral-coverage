# Helper functions
library(copula)
library(HAC)
library(methods)
library(truncnorm)
#library(tmvtnorm)

# ----------------------------------------
# Class stateAlgo
setClass("stateAlgo",
        slots = list(probaAcc="numeric",
                    lambda="numeric",
                    adaptMean="numeric",
                    adaptVar="numeric",
                    accRate="numeric"))


setGeneric (name= "adaptAlgo",
            def=function(object,iter,current_state,lambda,probaAcc,algo,prop,stochStep){standardGeneric("adaptAlgo")})

setMethod(f="adaptAlgo", 
          signature= c(object = "stateAlgo", iter = "numeric", current_state = "numeric", lambda = "numeric", probaAcc = "numeric", algo = "character", prop = "character"), 
          definition=function(object,iter,current_state,lambda,probaAcc,algo,prop,stochStep) 
          {
            x<-as.vector(current_state - object@adaptMean)
            object@adaptVar = object@adaptVar + (1/((iter)^stochStep)) * (diag(unlist(x)%*%t(unlist(x))) - object@adaptVar)
            object@adaptMean = object@adaptMean + (1/((iter)^stochStep))*unlist(x)
            
            # lambda is a scaling vector, of size 1 for MH algorithm, and of size P for Gibbs algorithms
            # the same is true for alpha, which is of size 1 for MH algos and of size P for Gibbs algos
            #print(paste("probaAcc in adaptAlgo is",probaAcc))
            object@probaAcc = probaAcc
            if (algo == "mh") {
              if(prop == "Gl") {
                object@lambda = object@lambda*exp((1/((iter)^stochStep))*(probaAcc - 0.234))
              } else if(prop == "Ad") {
                  object@lambda = 2.38*2.38/length(current_state)
                }
            } else if (prop == "CW") {
              for (j in 1:length(object@adaptMean))
              {
                if (probaAcc[j] > 0) # if alpha = 0 it means that this component was not randomly chosen by the HGS algo
                  object@lambda[j] = object@lambda[j]*exp((1/((iter)^stochStep))*(probaAcc[j] - 0.44))
              }
            } else if (prop == "Gl") {
              for (j in 1:length(object@adaptMean))
              {
                object@lambda[j]  = object@lambda[j]*exp((1/((iter)^stochStep))*(probaAcc[j] - 0.44));
              }
            }
            return(object)
          })
# ---------------------------------------

#' Params of Beta law from mode and sample size
#' @param mode the mode of the Beta dist
#' @param sample size the sample size
#' @return a and b, params of the Beta law
paramBetaFromModeSampSize<-function(m,s)
{
  a<-m*(s-2)+1
  b<-(s-2)*(1-m)+1
  return(list(a=a,b=b))
}

#' Params of Beta law from mode and sample size
#' @param mode the mode of the Beta dist
#' @param sample size the sample size
#' @return a and b, params of the Beta law
paramBetaFromMeanVar <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

#' Params of Beta law from similar triangular dist.
#' @param a min of the triangular dist.
#' @param b mode of the triangular dist.
#' @param c max of the triangular dist.
#' @return a and b, params of the Beta law
paramBetaFromTriang <- function(a,b,c) {
  m <- (a+b+c)/3
  v <- (a^2+b^2+c^2-a*c-a*b-b*c)/18
  res <- paramBetaFromMeanVar(m,v)
  return(params = res)
}

#' Density of the uniform proposal (only the part related to tau) computed as q(y|x)
#' @param y the point at which we want to evaluate the density
#' @param x the current point of the RW
#' @param sd the standard error 
jointDensTauUnif<-function(y,x,sd)
{
  binf<-max(0,x[1]-sqrt(3)*sd[1])
  bsup<-min(1,x[1]+sqrt(3)*sd[1])
  d1 <- dunif(y[1],min=binf,max=bsup)
  
  if(y[1] > x[2]){
    lambda<-min(sqrt(3)*sd[2],y[1]-x[2])
    binf<-max(0,x[2]-lambda)
    bsup<-min(1,x[2]+lambda)
    d2 <- dunif(y[2],min=binf,max=bsup)
  }else{
    d2 <- (y[1] == y[2]) 
  }
  
  return(d1*d2)
}

#' Density of the Gaussian proposal (only the part related to tau) computed as q(y|x)
#' @param y the point at which we want to evaluate the density
#' @param x the current point of the RW
#' @param sd the standard error 
jointDensTauNorm<-function(y,x,sd)
{
  d1 <- dtruncnorm(y[1],a=0,b=1,mean=x[1],sd=sd[1])
  
  if(y[1] > x[2]){
    d2 <- dtruncnorm(y[2],a=0,b=y[1],mean=x[2],sd=sd[2])
  }else{
    d2 <- (y[1] == y[2]) 
  }
  
  return(d1*d2)
}

#' Ratio of acceptance for random walk
#' Assume uniform priors
#' @param candidate the candidate generated by the MCMC procedure
#' @param current the current state of the chain
#' @param data the matrix of observations (assumed here to be of size n*4, and each row is modelled using a copula)
#' @param hyper the hyperparameters
#### !!! change hyper to more explicit notations?
ratioRW<-function(candidate,current,data,hyper,copFamily,dist,sdTau)
{
  #candidate<-unlist(candidate)
  #current<-unlist(current)
  size_state<-length(candidate)
  
  # ratio of priors
  ratioPriors = 1
  for (i in 1:(size_state-2))
  {
    ratioPriors = ratioPriors * (dunif(x=candidate[[i]],min=hyper[[i]][1],max=hyper[[i]][2]) / dunif(current[[i]],hyper[[i]][1],hyper[[i]][2]))
  }
  ratioPriors = ratioPriors * ((candidate["tauW"] >= 0) & (candidate["tauW"] <= 1)) 
  
  if (size_state > 3)
    ratioPriors =  ratioPriors * ((1/candidate["tauW"]) * ((candidate["tauB"] >= 0) & (candidate["tauB"] <= candidate["tauW"])))
  
  # ratio of likelihoods
  size_data = nrow(data)
  # we compute theta, the parameters of the copula, according to Kendall's tau values
  if (copFamily == "gumbel"){
    thetaWca<-copGumbel@iTau(candidate["tauW"])
    thetaWcu<-copGumbel@iTau(current$tauW)
    if (size_state > 3)
    {
      thetaBca<-copGumbel@iTau(candidate["tauB"])
      thetaBcu<-copGumbel@iTau(current$tauB)
    }
    typeHAC<-1
  }else if (copFamily == "clayton"){
    thetaWca<-copClayton@iTau(candidate["tauW"])
    thetaWcu<-copClayton@iTau(current$tauW)
    if (size_state > 3)
    {    
      thetaBca<-copClayton@iTau(candidate["tauB"])
      thetaBcu<-copClayton@iTau(current$tauB)  
    }
    typeHAC<-3
  }else if (copFamily == "frank"){
    thetaWca<-copFrank@iTau(candidate["tauW"])
    thetaWcu<-copFrank@iTau(current$tauW)
    if (size_state > 3)
    {
      thetaBca<-copFrank@iTau(candidate["tauB"])
      thetaBcu<-copFrank@iTau(current$tauB)    
    }
    typeHAC<-5
  }else if (copFamily == "amh"){
    thetaWca<-copAMH@iTau(candidate["tauW"])
    thetaWcu<-copAMH@iTau(current$tauW)
    if (size_state > 3)
    {
      thetaBca<-copAMH@iTau(candidate["tauB"])
      thetaBcu<-copAMH@iTau(current$tauB)    
    }
    typeHAC<-9
  }else if (copFamily == "joe"){
    thetaWca<-copJoe@iTau(candidate["tauW"])
    thetaWcu<-copJoe@iTau(current$tauW)
    if (size_state > 3)
    {
      thetaBca<-copJoe@iTau(candidate["tauB"])
      thetaBcu<-copJoe@iTau(current$tauB)    
    }
    typeHAC<-7
  }
  
  # define a and b param for beta distriutions according to mode and sample size
  # for the candidate
  a1_ca<-candidate["mu1"]*(candidate["s1"]-2)+1
  b1_ca<-(candidate["s1"]-2)*(1-candidate["mu1"])+1
  # for the current state
  a1_cu<-current$mu1*(current$s1-2)+1
  b1_cu<-(current$s1-2)*(1-current$mu1)+1
  
  # we define the two copulas, one for the candidate param and one for the current param
  mat<-as.matrix(data)
  if (size_state > 3)
  {
    copCand<-hac(type=typeHAC,tree=list(list("y1","y2",thetaWca),list("y3","y4",thetaWca),thetaBca))
    copCurr<-hac(type=typeHAC,tree=list(list("y1","y2",thetaWcu),list("y3","y4",thetaWcu),thetaBcu))
    
    a2_ca<-candidate["mu2"]*(candidate["s2"]-2)+1
    b2_ca<-(candidate["s2"]-2)*(1-candidate["mu2"])+1
    a2_cu<-current$mu2*(current$s2-2)+1
    b2_cu<-(current$s2-2)*(1-current$mu2)+1
    
    marginMat_ca <- mat * 0
    marginMat_ca[,1] <- pbeta(mat[,1],a1_ca,b1_ca)
    marginMat_ca[,2] <- pbeta(mat[,2],a1_ca,b1_ca)
    marginMat_ca[,3] <- pbeta(mat[,3],a2_ca,b2_ca)
    marginMat_ca[,4] <- pbeta(mat[,4],a2_ca,b2_ca)    
    marginMat_cu <- mat * 0
    marginMat_cu[,1] <- pbeta(mat[,1],a1_ca,b1_cu)
    marginMat_cu[,2] <- pbeta(mat[,2],a1_ca,b1_cu)
    marginMat_cu[,3] <- pbeta(mat[,3],a2_ca,b2_cu)
    marginMat_cu[,4] <- pbeta(mat[,4],a2_ca,b2_cu)    
    
    num<-dHAC(marginMat_ca,copCand)*dbeta(mat[,1],a1_ca,b1_ca)*dbeta(mat[,2],a1_ca,b1_ca)*dbeta(mat[,3],a2_ca,b2_ca)*dbeta(mat[,4],a2_ca,b2_ca)
    denom<-dHAC(marginMat_cu,copCurr)*dbeta(mat[,1],a1_cu,b1_cu)*dbeta(mat[,2],a1_cu,b1_cu)*dbeta(mat[,3],a2_cu,b2_cu)*dbeta(mat[,4],a2_cu,b2_cu)
    
  } else {
    copCand<-hac(type=typeHAC,tree=list("y1","y2",thetaWca))
    copCurr<-hac(type=typeHAC,tree=list("y1","y2",thetaWcu))
    
    marginMat_ca <- mat * 0
    marginMat_ca[,1] <- pbeta(mat[,1],a1_ca,b1_ca)
    marginMat_ca[,2] <- pbeta(mat[,2],a1_ca,b1_ca)    
    marginMat_cu <- mat * 0
    marginMat_cu[,1] <- pbeta(mat[,1],a1_cu,b1_cu)
    marginMat_cu[,2] <- pbeta(mat[,2],a1_cu,b1_cu)      
    
    num<-dHAC(marginMat_ca,copCand)*dbeta(mat[,1],a1_ca,b1_ca)*dbeta(mat[,2],a1_ca,b1_ca)
    denom<-dHAC(marginMat_cu,copCurr)*dbeta(mat[,1],a1_cu,b1_cu)*dbeta(mat[,2],a1_cu,b1_cu)
  }
  
  #print(prod(na.omit(denom)))
  ratioLik = exp(sum(log(num)) - sum(log(denom)))  #prod(na.omit(num))/prod(na.omit(denom))
  if (prod(na.omit(num)) == 0 & prod(na.omit(denom)) == 0)
    ratioLik = 0
  if (any(is.na(num)))
  {
    print(num)
    ratioLik = 0
  }
  
  # ratio of proposals 
  if (size_state > 3)
  {
    cu<-c(current$tauW,current$tauB)
    ca<-c(candidate["tauW"],candidate["tauB"])
    if (dist == "unif"){
      ratioProp = jointDensTauUnif(cu,ca,sdTau) / jointDensTauUnif(ca,cu,sdTau)
    }else if (dist == "norm"){
      ratioProp = (jointDensTauNorm(cu,ca,sdTau) / jointDensTauNorm(ca,cu,sdTau)) * (dtruncnorm(current$mu1, a = 0, b = 1) / dtruncnorm(candidate["mu1"], a = 0, b = 1)) * (dtruncnorm(current$mu2, a = 0, b = 1) / dtruncnorm(candidate["mu2"], a = 0, b = 1))
    }
    if (is.nan(ratioProp[[1]]))
      ratioProp<-0
  }
  else {ratioProp <- 1}
  
  #print(paste("ratio priors :",ratioPriors,"ratio likelihood :",ratioLik,"ratioProp :",ratioProp))
  return(ratioPriors*ratioLik*ratioProp)
}


#' Ratio of acceptance for OSR (only 1 period)
#' Assume uniform priors
#' @param candidate the candidate generated by the MCMC procedure
#' @param current the current state of the chain
#' @param data the matrix of observations (assumed here to be of size n*4, and each row is modelled using a copula)
#' @param hyper the hyperparameters
#### !!! change hyper to more explicit notations?
ratioRW_osr<-function(candidate,current,data,hyper,copFamily,dist,sdTau)
{
  #candidate<-unlist(candidate)
  #current<-unlist(current)
  size_state<-length(candidate)
  
  # ratio of priors
  ratioPriors = 1
  for (i in 1:(size_state-1))
  {
    ratioPriors = ratioPriors * (dunif(x=candidate[[i]],min=hyper[[i]][1],max=hyper[[i]][2]) / dunif(current[[i]],hyper[[i]][1],hyper[[i]][2]))
  }
  ratioPriors = ratioPriors * ((candidate["tauW"] >= 0) & (candidate["tauW"] <= 1)) 
  
  # ratio of likelihoods
  size_data = nrow(data)
  # we compute theta, the parameters of the copula, according to Kendall's tau values
  if (copFamily == "gumbel"){
    thetaWca<-copGumbel@iTau(candidate["tauW"])
    thetaWcu<-copGumbel@iTau(current$tauW)
    typeHAC<-1
  }else if (copFamily == "clayton"){
    thetaWca<-copClayton@iTau(candidate["tauW"])
    thetaWcu<-copClayton@iTau(current$tauW)
    typeHAC<-3
  }else if (copFamily == "frank"){
    thetaWca<-copFrank@iTau(candidate["tauW"])
    thetaWcu<-copFrank@iTau(current$tauW)
    typeHAC<-5
  }else if (copFamily == "amh"){
    thetaWca<-copAMH@iTau(candidate["tauW"])
    thetaWcu<-copAMH@iTau(current$tauW)
    typeHAC<-9
  }else if (copFamily == "joe"){
    thetaWca<-copJoe@iTau(candidate["tauW"])
    thetaWcu<-copJoe@iTau(current$tauW)
    typeHAC<-7
  }
  
  # define a and b param for beta distriutions according to mode and sample size
  # for the candidate
  a1_ca<-candidate["mu1"]*(candidate["s1"]-2)+1
  b1_ca<-(candidate["s1"]-2)*(1-candidate["mu1"])+1
  # for the current state
  a1_cu<-current$mu1*(current$s1-2)+1
  b1_cu<-(current$s1-2)*(1-current$mu1)+1
  
  # we define the two copulas, one for the candidate param and one for the current param
  mat<-as.matrix(data)
  copCand<-hac(type=typeHAC,tree=list("y1","y2",thetaWca))
  copCurr<-hac(type=typeHAC,tree=list("y1","y2",thetaWcu))

  marginMat_ca <- mat * 0
  marginMat_ca[,1] <- pbeta(mat[,1],a1_ca,b1_ca)
  marginMat_ca[,2] <- pbeta(mat[,2],a1_ca,b1_ca)
  marginMat_cu <- mat * 0
  marginMat_cu[,1] <- pbeta(mat[,1],a1_cu,b1_cu)
  marginMat_cu[,2] <- pbeta(mat[,2],a1_cu,b1_cu)
  
  num<-dHAC(marginMat_ca,copCand)*dbeta(mat[,1],a1_ca,b1_ca)*dbeta(mat[,2],a1_ca,b1_ca)
  denom<-dHAC(marginMat_cu,copCurr)*dbeta(mat[,1],a1_cu,b1_cu)*dbeta(mat[,2],a1_cu,b1_cu)
  
  ratioLik = prod(na.omit(num))/prod(na.omit(denom))
  if (prod(na.omit(num)) == 0 & prod(na.omit(denom)) == 0)
    ratioLik = 0
  
  # ratio of proposals (only for normal proposal since in the uniform case the densities are equal to 1)
  ratioProp = 1
  if (dist == "norm"){
  ratioProp = ratioProp * dtruncnorm(candidate["mu1"], a = 0, b = 1) * dtruncnorm(candidate["tauW"], a = 0, b = 1) / 
              (dtruncnorm(current$mu1, a = 0, b = 1) * dtruncnorm(current$tauW, a = 0, b = 1))
  }
  
  #print(paste("ratio priors :",ratioPriors,"ratio likelihood :",ratioLik,"ratioProp :",ratioProp))
  return(ratioPriors*ratioLik*ratioProp)
}

#' Ratio of acceptance for independence case
#' Assume uniform priors
#' @param candidate the candidate generated by the MCMC procedure
#' @param current the current state of the chain
#' @param data the matrix of observations (assumed here to be of size n*4, and each row is modelled using a copula)
#' @param hyper the hyperparameters
#### !!! change hyper to more explicit notations?
ratioRW_indep<-function(candidate,current,data,hyper,dist)
{
  #candidate<-unlist(candidate)
  #current<-unlist(current)
  size_state<-length(candidate)
  
  # ratio of priors
  ratioPriors = 1
  for (i in 1:size_state)
  {
    ratioPriors = ratioPriors * (dunif(x=candidate[[i]],min=hyper[[i]][1],max=hyper[[i]][2]) / dunif(current[[i]],hyper[[i]][1],hyper[[i]][2]))
  }
  
  # ratio of likelihoods
  size_data = nrow(data)

  # define a and b param for beta distriutions according to mode and sample size
  # for the candidate
  a1_ca<-candidate["mu1"]*(candidate["s1"]-2)+1
  b1_ca<-(candidate["s1"]-2)*(1-candidate["mu1"])+1
  a2_ca<-candidate["mu2"]*(candidate["s2"]-2)+1
  b2_ca<-(candidate["s2"]-2)*(1-candidate["mu2"])+1
  # for the current state
  a1_cu<-current$mu1*(current$s1-2)+1
  b1_cu<-(current$s1-2)*(1-current$mu1)+1
  a2_cu<-current$mu2*(current$s2-2)+1
  b2_cu<-(current$s2-2)*(1-current$mu2)+1
  
  num<-dbeta(mat[,1],a1_ca,b1_ca)*dbeta(mat[,2],a1_ca,b1_ca)*dbeta(mat[,3],a2_ca,b2_ca)*dbeta(mat[,4],a2_ca,b2_ca)
  denom<-dbeta(mat[,1],a1_cu,b1_cu)*dbeta(mat[,2],a1_cu,b1_cu)*dbeta(mat[,3],a2_cu,b2_cu)*dbeta(mat[,4],a2_cu,b2_cu)
  
  ratioLik = prod(na.omit(num))/prod(na.omit(denom))
  if (prod(na.omit(num)) == 0 & prod(na.omit(denom)) == 0)
    ratioLik = 0
  
  # ratio of proposals
  ratioProp = 1;
  if (dist == "norm"){
    ratioProp = ratioProp * dtruncnorm(candidate["mu1"], a = 0, b = 1) / dtruncnorm(current$mu1, a = 0, b = 1)
  }
  
  #print(paste("ratio priors :",ratioPriors,"ratio likelihood :",ratioLik,"ratioProp :",ratioProp))
  return(ratioPriors*ratioLik*ratioProp)
}

#' Ratio of acceptance for independence case and OSR
#' Assume uniform priors
#' @param candidate the candidate generated by the MCMC procedure
#' @param current the current state of the chain
#' @param data the matrix of observations (assumed here to be of size n*4, and each row is modelled using a copula)
#' @param hyper the hyperparameters
#### !!! change hyper to more explicit notations?
ratioRW_osr_indep<-function(candidate,current,data,hyper)
{
  #candidate<-unlist(candidate)
  #current<-unlist(current)
  size_state<-length(candidate)
  
  # ratio of priors
  ratioPriors = 1
  for (i in 1:size_state)
  {
    ratioPriors = ratioPriors * (dunif(x=candidate[[i]],min=hyper[[i]][1],max=hyper[[i]][2]) / dunif(current[[i]],hyper[[i]][1],hyper[[i]][2]))
  }
  
  # ratio of likelihoods
  size_data = nrow(data)
  
  # define a and b param for beta distriutions according to mode and sample size
  # for the candidate
  a1_ca<-candidate["mu1"]*(candidate["s1"]-2)+1
  b1_ca<-(candidate["s1"]-2)*(1-candidate["mu1"])+1
  # for the current state
  a1_cu<-current$mu1*(current$s1-2)+1
  b1_cu<-(current$s1-2)*(1-current$mu1)+1
  
  num<-dbeta(mat[,1],a1_ca,b1_ca)*dbeta(mat[,2],a1_ca,b1_ca)
  denom<-dbeta(mat[,1],a1_cu,b1_cu)*dbeta(mat[,2],a1_cu,b1_cu)
  
  ratioLik = exp(sum(log(num)) - sum(log(denom)))  #prod(na.omit(num))/prod(na.omit(denom))
  if (prod(na.omit(num)) == 0 & prod(na.omit(denom)) == 0)
    ratioLik = 0
  if (any(is.na(num)))
  {
    print(num)
    ratioLik = 0
  }
  
  # ratio of proposals
  ratioProp = 1;
  if (dist == "norm"){
    ratioProp = ratioProp * dtruncnorm(candidate["mu1"], a = 0, b = 1) / dtruncnorm(current$mu1, a = 0, b = 1)
  }
  
  #print(paste("ratio priors :",ratioPriors,"ratio likelihood :",ratioLik,"ratioProp :",ratioProp))
  return(ratioPriors*ratioLik*ratioProp)
}


#' Uniform proposal for the last two components as they are not independent from each other
#' @param last the last state of the chain (used to compute the random walk on the first component)
#' @param sd the standard error for the random walk
jointPropTauUnif<-function(last,sd)
{
  # Generate the first component according to p(x)
  # Truncate the proposal between 0 and 1
  binf<-max(0,last[1]-sqrt(3)*sd[1])
  bsup<-min(1,last[1]+sqrt(3)*sd[1])
  u1<-runif(1,min=binf,max=bsup)
  # Generate the second component accoding to p(y|x)
  if (last[2] < u1){
    lambda<-min(sqrt(3)*sd[2],u1-last[2])
    binf<-max(0,last[2]-lambda)
    bsup<-min(1,last[2]+lambda)
    u2<-runif(1,min=binf,max=bsup)
  }else{
    u2<-u1
  }
  return(c(u1,u2))
}

#' Gaussian proposal for the last two components as they are not independent from each other
#' @param last the last state of the chain (used to compute the random walk on the first component)
#' @param sd the standard error for the random walk
jointPropTauNorm<-function(last,sd)
{
  # Generate the first component according to p(x)
  # Truncate the proposal between 0 and 1 (use the truncnorm package)
  u1<-rtruncnorm(1,a=0,b=1,mean=last[1],sd=sd[1])
  # Generate the second component according to p(y|x)
  if (last[2] < u1){
    u2<-rtruncnorm(1,a=0,b=u1,mean=last[2],sd=sd[2])
  }else{
    u2<-u1
  }
  return(c(u1,u2))
}


#' Update the acceptance rate of the chain 
#' @param chain the chain
#' @param iter the current iteration
#' @param lastAR the last value of the aceptance rate
updateAccRate<-function(chain,iter,lastAR) 
{
  hasMoved<-(chain[iter,] != chain[(iter-1),])
  accRate <- (lastAR*(iter-1) + hasMoved)/iter
  return(accRate)
}


#' Generate a MCMC sample
#' @param last_state the last state of the chain
#' @param state_algo an object of type stateAlgo
#' @param algo the type of algorithm ("mh" or "hgs")
#' @param prop the proposal
#' @param hyper the vector of hyperparameters
#' @param data the matrix of data
generateMCMC<-function(last_state,state_algo,algo,prop,dist,hyper,copFamily,data)
{
  size_state<-length(state_algo@adaptMean)
  
  next_state<-rep(0,size_state)
  candidate<-rep(0,size_state)
  names(candidate)<-c("mu1","mu2","s1","s2","tauW","tauB")
  
  sd<-sqrt(state_algo@lambda*state_algo@adaptVar)   

  if(algo == "mh")
  {
    a<-last_state-sqrt(3)*sd
    b<-last_state+sqrt(3)*sd
    a[[1]]<-max(0,a[[1]]);a[[2]]<-max(0,a[[2]]); # troncate mode to be between 0 and 1
    b[[1]]<-min(1,b[[1]]);b[[2]]<-min(1,b[[2]]); # troncate mode to be between 0 and 1
    for (j in 1:(size_state-2)) # generate candidate for all component except the last one which needs to be computed differently (see below)
    {
      if (dist == "unif"){
        # compute the bounds of the uniform random walk given the standard error and the current state -> we want a uniform distribution
        # centered on the current state and of a given variance, computed with the adaptive scheme
        candidate[j] <- runif(1,min=a[[j]],max=b[[j]])
      }else if(dist == "norm"){
        if(j <= 2){
          candidate[j] <- rtruncnorm(1,a=0,b=1,mean=last_state[[j]],sd=sd[j])
        }else{
         candidate[j] <- rnorm(1,mean=last_state[[j]],sd=sd[j])
        }
      }
    }
    # the last components are the copula parameters and they should be sampled from jointly
    tauLast<-c(last_state$tauW,last_state$tauB)
    if (dist == "unif"){
      tau<-jointPropTauUnif(tauLast,sd[(size_state-1):size_state])
      candidate[(size_state-1):size_state]<-tau
    }else if(dist == "norm"){
      tau<-jointPropTauNorm(tauLast,sd[(size_state-1):size_state])
      candidate[(size_state-1):size_state]<-tau
    }
    
    #print(candidate)
    r = min(1,ratioRW(candidate,last_state,data,hyper,copFamily,dist,sd[(size_state-1):size_state]))
    probaAcc<-r
    #print(probaAcc)
    
    u<-runif(1)
    if(is.finite(r) & u <= r) {
      next_state<-candidate
    } else {
      next_state<-last_state
    }
  }
  else if(algo == "hgs")
  {
    if(prop == "CW")
    {
      candidate = unlist(last_state[1,])
      # we randomly choose one of the components  
      p = sample(1:(size_state-1),1)
      if(p <= 2){
        if (dist == "unif"){
          a<-max(0,unlist(last_state[p]-sqrt(3)*sd[p]))
          b<-min(1,unlist(last_state[p]+sqrt(3)*sd[p]))
          candidate[p] <- runif(1,min=a,max=b)
        }else if(dist == "norm"){
          candidate[p] <- rtruncnorm(1,a=0,b=1,mean=last_state[[p]],sd=sd[p])
        }
      }else if(p < (size_state-1)){
        if (dist == "unif"){
          a<-last_state[p]-sqrt(3)*sd[p]
          b<-last_state[p]+sqrt(3)*sd[p]
          candidate[p] <- runif(1,min=a[[1]],max=b[[1]])
        }else if(dist == "norm"){
          candidate[p] <- rnorm(1,mean=last_state[[p]],sd=sd[p])
        }
      }else{
        tauLast<-c(last_state$tauW,last_state$tauB)
        if (dist == "unif"){
          tau<-jointPropTauUnif(tauLast,sd[(size_state-1):size_state])
          candidate[(size_state-1):size_state]<-tau
        }else if(dist == "norm"){
          tau<-jointPropTauNorm(tauLast,sd[(size_state-1):size_state])
          candidate[(size_state-1):size_state]<-tau
        }
      }
      #print(paste("candidate :",candidate))
    }
    else if(prop == "Gl") ## check the differences with mh sampling ... mh should be multivariate?
    {
      a<-last_state-sqrt(3)*sd
      b<-last_state+sqrt(3)*sd
      a[[1]]<-max(0,a[[1]]);a[[2]]<-max(0,a[[2]]); # troncate mode to be between 0 and 1
      b[[1]]<-min(1,b[[1]]);b[[2]]<-min(1,b[[2]]); # troncate mode to be between 0 and 1      
      for (p in 1:(size_state-1))
      {
        if (dist == "unif"){
          candidate[p] <- runif(1,min=a[[p]],max=b[[p]])
        }else if (dist == "norm")
          if(p <= 2){
            candidate[p] <- rtruncnorm(1,a=0,b=1,mean=last_state[[p]],sd=sd[[p]])
          }else{
            candidate[p] <- rnorm(1,mean=last_state[[p]],sd=sd[[p]])
        }
      }
      tauLast<-c(last_state$tauW,last_state$tauB)
      if (dist == "unif"){
        tau<-jointPropTauUnif(tauLast,sd[(size_state-1):size_state])
        candidate[(size_state-1):size_state]<-tau
      }else if(dist == "norm"){
        tau<-jointPropTauNorm(tauLast,sd[(size_state-1):size_state])
        candidate[(size_state-1):size_state]<-tau
      }
    }

    r = min(1,ratioRW(candidate,last_state,data,hyper,copFamily,dist,sd[(size_state-1):size_state]))
    
    u<-runif(1)
    if(is.finite(r) & u <= r) {
      next_state<-candidate
    }else{
      next_state<-last_state
    }
    
    probaAcc<-rep(0,size_state)  
    if(prop == "CW")
    {
      probaAcc[p]<-r  
    }
    else if(prop == "Gl")
    {
      temp<-last_state
      for (j in 1:size_state)
      {
        temp[j] = candidate[j]
        r = min(1,ratioRW(candidate,last_state,data,hyper,copFamily,dist,sd[(size_state-1):size_state]))
        
        probaAcc[j] = r
        temp[j] = last_state[j]
      }
    }
  }
  #print(paste("r =",r,"probaAcc =",probaAcc))
  return(list(chain=next_state,lambda=lambda,accRate=probaAcc))
}


#' Generate MCMC for OSR data (less observed dates so the code must be adapted and to enhance clarity and avoid multiple if statements we use
#' a different function)
generateMCMC_OSR<-function(last_state,state_algo,algo,prop,dist,hyper,copFamily,data)
{
  size_state<-length(state_algo@adaptMean)
  
  next_state<-rep(0,size_state)
  candidate<-rep(0,size_state)
  names(candidate)<-c("mu1","s1","tauW")
  
  sd<-sqrt(state_algo@lambda*state_algo@adaptVar)   
  
  if(algo == "mh")
  {
    a<-last_state-sqrt(3)*sd
    b<-last_state+sqrt(3)*sd
    a[[1]]<-max(0,a[[1]]);a[[3]]<-max(0,a[[3]]); # troncate mode to be between 0 and 1
    b[[1]]<-min(1,b[[1]]);b[[3]]<-min(1,b[[3]]); # troncate tau to be between 0 and 1
    for (j in 1:size_state) 
    {
      if (dist == "unif"){
        # compute the bounds of the uniform random walk given the standard error and the current state -> we want a uniform distribution
        # centered on the current state and of a given variance, computed with the adaptive scheme
        candidate[j] <- runif(1,min=a[[j]],max=b[[j]])
      }else if(dist == "norm"){
        if(j != 2){
          candidate[j] <- rtruncnorm(1,a=0,b=1,mean=last_state[[j]],sd=sd[j])
        }else{
          candidate[j] <- rnorm(1,mean=last_state[[j]],sd=sd[j])
        }
      }
    }
    
    #print(candidate)
    r = min(1,ratioRW_osr(candidate,last_state,data,hyper,copFamily,dist,sd[(size_state-1):size_state]))
    probaAcc<-r
    #print(probaAcc)
    
    u<-runif(1)
    if(is.finite(r) & u <= r) {
      next_state<-candidate
    } else {
      next_state<-last_state
    }
  }
  else if(algo == "hgs")
  {
    if(prop == "CW")
    {
      candidate = unlist(last_state[1,])
      # we randomly choose one of the components  
      p = sample(1:size_state,1)
      if(p != 2){
        if (dist == "unif"){
          a<-max(0,unlist(last_state[p]-sqrt(3)*sd[p]))
          b<-min(1,unlist(last_state[p]+sqrt(3)*sd[p]))
          candidate[p] <- runif(1,min=a,max=b)
        }else if(dist == "norm"){
          candidate[p] <- rtruncnorm(1,a=0,b=1,mean=last_state[[p]],sd=sd[p])
        }
      }else{
        if (dist == "unif"){
          a<-last_state[p]-sqrt(3)*sd[p]
          b<-last_state[p]+sqrt(3)*sd[p]
          candidate[p] <- runif(1,min=a[[1]],max=b[[1]])
        }else if(dist == "norm"){
          candidate[p] <- rnorm(1,mean=last_state[[p]],sd=sd[p])
        }
      }
      #print(paste("candidate :",candidate))
    }
    else if(prop == "Gl") ## check the differences with mh sampling ... mh should be multivariate?
    {
      a<-last_state-sqrt(3)*sd
      b<-last_state+sqrt(3)*sd
      a[[1]]<-max(0,a[[1]]);a[[3]]<-max(0,a[[3]]); # troncate mode to be between 0 and 1
      b[[1]]<-min(1,b[[1]]);b[[3]]<-min(1,b[[3]]); # troncate tau to be between 0 and 1   
      for (p in 1:size_state)
      {
        if (dist == "unif"){
          candidate[p] <- runif(1,min=a[[p]],max=b[[p]])
        }else if (dist == "norm")
          if(p != 2){
            candidate[p] <- rtruncnorm(1,a=0,b=1,mean=last_state[[p]],sd=sd[[p]])
          }else{
            candidate[p] <- rnorm(1,mean=last_state[[p]],sd=sd[[p]])
          }
      }
    }
    
    r = min(1,ratioRW_osr(candidate,last_state,data,hyper,copFamily,dist,sd[(size_state-1):size_state]))
    
    u<-runif(1)
    if(is.finite(r) & u <= r) {
      next_state<-candidate
    }else{
      next_state<-last_state
    }
    
    probaAcc<-rep(0,size_state)  
    if(prop == "CW")
    {
      probaAcc[p]<-r  
    }
    else if(prop == "Gl")
    {
      temp<-last_state
      for (j in 1:size_state)
      {
        temp[j] = candidate[j]
        r = min(1,ratioRW_osr(candidate,last_state,data,hyper,copFamily,dist,sd[(size_state-1):size_state]))
        
        probaAcc[j] = r
        temp[j] = last_state[j]
      }
    }
  }
  #print(paste("r =",r,"probaAcc =",probaAcc))
  return(list(chain=next_state,lambda=lambda,accRate=probaAcc))
}


#' Generate a MCMC sample in the independent case
#' @param last_state the last state of the chain
#' @param state_algo an object of type stateAlgo
#' @param algo the type of algorithm ("mh" or "hgs")
#' @param prop the proposal
#' @param hyper the vector of hyperparameters
#' @param data the matrix of data
generateMCMC_indep<-function(last_state,state_algo,algo,prop,dist,hyper,data)
{
  size_state<-length(state_algo@adaptMean)
  
  next_state<-rep(0,size_state)
  candidate<-rep(0,size_state)
  names(candidate)<-c("mu1","mu2","s1","s2")
  
  sd<-sqrt(state_algo@lambda*state_algo@adaptVar)   
  
  if(algo == "mh")
  {
    a<-last_state-sqrt(3)*sd
    b<-last_state+sqrt(3)*sd
    a[[1]]<-max(0,a[[1]]);a[[2]]<-max(0,a[[2]]); # troncate mode to be between 0 and 1
    b[[1]]<-min(1,b[[1]]);b[[2]]<-min(1,b[[2]]); # troncate mode to be between 0 and 1
    for (j in 1:size_state) # generate candidate for all component except the last one which needs to be computed differently (see below)
    {
      if (dist == "unif"){
        # compute the bounds of the uniform random walk given the standard error and the current state -> we want a uniform distribution
        # centered on the current state and of a given variance, computed with the adaptive scheme
        candidate[j] <- runif(1,min=a[[j]],max=b[[j]])
      }else if(dist == "norm"){
        if(j <= 2){
          candidate[j] <- rtruncnorm(1,a=0,b=1,mean=last_state[[j]],sd=sd[j])
        }else{
          candidate[j] <- rnorm(1,mean=last_state[[j]],sd=sd[j])
        }
      }
    }

    #print(candidate)
    r = min(1,ratioRW_indep(candidate,last_state,data,hyper,dist))
    probaAcc<-r
    #print(probaAcc)
    
    u<-runif(1)
    if(is.finite(r) & u <= r) {
      next_state<-candidate
    } else {
      next_state<-last_state
    }
  }
  else if(algo == "hgs")
  {
    if(prop == "CW")
    {
      candidate = unlist(last_state[1,])
      # we randomly choose one of the components  
      p = sample(1:size_state,1)
      if(p <= 2){
        if (dist == "unif"){
          a<-max(0,unlist(last_state[p]-sqrt(3)*sd[p]))
          b<-min(1,unlist(last_state[p]+sqrt(3)*sd[p]))
          candidate[p] <- runif(1,min=a,max=b)
        }else if(dist == "norm"){
          candidate[p] <- rtruncnorm(1,a=0,b=1,mean=last_state[[p]],sd=sd[p])
        }
      }else{
        if (dist == "unif"){
          a<-last_state[p]-sqrt(3)*sd[p]
          b<-last_state[p]+sqrt(3)*sd[p]
          candidate[p] <- runif(1,min=a[[1]],max=b[[1]])
        }else if(dist == "norm"){
          candidate[p] <- rnorm(1,mean=last_state[[p]],sd=sd[p])
        }
      }
      #print(paste("candidate :",candidate))
    }
    else if(prop == "Gl") ## check the differences with mh sampling ... mh should be multivariate?
    {
      a<-last_state-sqrt(3)*sd
      b<-last_state+sqrt(3)*sd
      a[[1]]<-max(0,a[[1]]);a[[2]]<-max(0,a[[2]]); # troncate mode to be between 0 and 1
      b[[1]]<-min(1,b[[1]]);b[[2]]<-min(1,b[[2]]); # troncate mode to be between 0 and 1      
      for (p in 1:size_state)
      {
        if (dist == "unif"){
          candidate[p] <- runif(1,min=a[[p]],max=b[[p]])
        }else if (dist == "norm")
          if(p <= 2){
            candidate[p] <- rtruncnorm(1,a=0,b=1,mean=last_state[[p]],sd=sd[[p]])
          }else{
            candidate[p] <- rnorm(1,mean=last_state[[p]],sd=sd[[p]])
          }
      }
    }
    
    r = min(1,ratioRW_indep(candidate,last_state,data,hyper,dist))
    
    u<-runif(1)
    if(is.finite(r) & u <= r) {
      next_state<-candidate
    }else{
      next_state<-last_state
    }
    
    probaAcc<-rep(0,size_state)  
    if(prop == "CW")
    {
      probaAcc[p]<-r  
    }
    else if(prop == "Gl")
    {
      temp<-last_state
      for (j in 1:size_state)
      {
        temp[j] = candidate[j]
        r = min(1,ratioRW_indep(candidate,last_state,data,hyper,dist))
        
        probaAcc[j] = r
        temp[j] = last_state[j]
      }
    }
  }
  #print(paste("r =",r,"probaAcc =",probaAcc))
  return(list(chain=next_state,lambda=lambda,accRate=probaAcc))
}


#' Generate a MCMC sample in the independent case and for OSR
#' @param last_state the last state of the chain
#' @param state_algo an object of type stateAlgo
#' @param algo the type of algorithm ("mh" or "hgs")
#' @param prop the proposal
#' @param hyper the vector of hyperparameters
#' @param data the matrix of data
generateMCMC_osr_indep<-function(last_state,state_algo,algo,prop,dist,hyper,data)
{
  size_state<-length(state_algo@adaptMean)
  
  next_state<-rep(0,size_state)
  candidate<-rep(0,size_state)
  names(candidate)<-c("mu1","s1")
  
  sd<-sqrt(state_algo@lambda*state_algo@adaptVar)   
  
  if(algo == "mh")
  {
    a<-last_state-sqrt(3)*sd
    b<-last_state+sqrt(3)*sd
    a[[1]]<-max(0,a[[1]])# troncate mode to be between 0 and 1
    b[[1]]<-min(1,b[[1]])
    for (j in 1:size_state) 
    {
      if (dist == "unif"){
        # compute the bounds of the uniform random walk given the standard error and the current state -> we want a uniform distribution
        # centered on the current state and of a given variance, computed with the adaptive scheme
        candidate[j] <- runif(1,min=a[[j]],max=b[[j]])
      }else if(dist == "norm"){
        if(j == 1){
          candidate[j] <- rtruncnorm(1,a=0,b=1,mean=last_state[[j]],sd=sd[j])
        }else{
          candidate[j] <- rnorm(1,mean=last_state[[j]],sd=sd[j])
        }
      }
    }
    
    #print(candidate)
    r = min(1,ratioRW_osr_indep(candidate,last_state,data,hyper))
    probaAcc<-r
    #print(probaAcc)
    
    u<-runif(1)
    if(is.finite(r) & u <= r) {
      next_state<-candidate
    } else {
      next_state<-last_state
    }
  }
  else if(algo == "hgs")
  {
    if(prop == "CW")
    {
      candidate = unlist(last_state[1,])
      # we randomly choose one of the components  
      p = sample(1:size_state,1)
      if(p == 1){
        if (dist == "unif"){
          a<-max(0,unlist(last_state[p]-sqrt(3)*sd[p]))
          b<-min(1,unlist(last_state[p]+sqrt(3)*sd[p]))
          candidate[p] <- runif(1,min=a,max=b)
        }else if(dist == "norm"){
          candidate[p] <- rtruncnorm(1,a=0,b=1,mean=last_state[[p]],sd=sd[p])
        }
      }else{
        if (dist == "unif"){
          a<-last_state[p]-sqrt(3)*sd[p]
          b<-last_state[p]+sqrt(3)*sd[p]
          candidate[p] <- runif(1,min=a[[1]],max=b[[1]])
        }else if(dist == "norm"){
          candidate[p] <- rnorm(1,mean=last_state[[p]],sd=sd[p])
        }
      }
      #print(paste("candidate :",candidate))
    }
    else if(prop == "Gl") ## check the differences with mh sampling ... mh should be multivariate?
    {
      a<-last_state-sqrt(3)*sd
      b<-last_state+sqrt(3)*sd
      a[[1]]<-max(0,a[[1]]); # troncate mode to be between 0 and 1
      b[[1]]<-min(1,b[[1]]); # troncate tau to be between 0 and 1   
      for (p in 1:size_state)
      {
        if (dist == "unif"){
          candidate[p] <- runif(1,min=a[[p]],max=b[[p]])
        }else if (dist == "norm")
          if(p == 1){
            candidate[p] <- rtruncnorm(1,a=0,b=1,mean=last_state[[p]],sd=sd[[p]])
          }else{
            candidate[p] <- rnorm(1,mean=last_state[[p]],sd=sd[[p]])
          }
      }
    }
    
    r = min(1,ratioRW_osr_indep(candidate,last_state,data,hyper))
    
    u<-runif(1)
    if(is.finite(r) & u <= r) {
      next_state<-candidate
    }else{
      next_state<-last_state
    }
    
    probaAcc<-rep(0,size_state)  
    if(prop == "CW")
    {
      probaAcc[p]<-r  
    }
    else if(prop == "Gl")
    {
      temp<-last_state
      for (j in 1:size_state)
      {
        temp[j] = candidate[j]
        r = min(1,ratioRW_osr_indep(candidate,last_state,data,hyper))
        
        probaAcc[j] = r
        temp[j] = last_state[j]
      }
    }
  }
  #print(paste("r =",r,"probaAcc =",probaAcc))
  return(list(chain=next_state,lambda=lambda,accRate=probaAcc))
}

#' Batch means computation
#' @param chain the Markov Chain univariate
#' @return the estimated variance of the different components of the chain
batchVar<-function(chain)
{
  n<-length(chain)
  batchSize<-floor(sqrt(n))
  
  nbBatchs<-(n-batchSize+1)
  meanBatchs<-rep(0,nbBatchs)
  meanChain<-mean(chain)
  for (i in 1:nbBatchs)
  {
    meanBatchs[i]<-mean(chain[i:(i+batchSize-1)])
  }
  varChain<-sum((meanBatchs-meanChain)^2) * (n*batchSize/((n-batchSize)*(n-batchSize+1)))
  return(varChain)
}


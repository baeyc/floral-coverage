# Floral coverage model

#options(echo=FALSE)
#args<-commandArgs(trailingOnly = TRUE)
#print(args)

args <- c("pasture","4","0","0.7","0","0.7","10","1")

landuseCat <-args[1]
nbrep <- as.integer(args[2])
m1l <- as.numeric(args[3])
m1u <- as.numeric(args[4])
m2l <- as.numeric(args[5])
m2u <- as.numeric(args[6])
nb_mcmc <- as.numeric(args[7])
nb_burn <- as.numeric(args[8])
rm(args)

setwd("/home/baeyc/Code/es2cropBB/pollination/calibration/floralCoverage/")

# load libraries
library(R.oo)
library(gridExtra)
library(reshape)

source("helpers.R")
source("multiplots.R")

# load data -> we have one matrix of observations per landscape category
source("loadObsData.R")
data<-eval(parse(text=landuseCat))

data <- na.omit(data)

## Priors
# 1. Hyperparameters (to be changed for each landscape category)
# prior on the mode and sample size of the beta law for the first period
hypMu1<-c(m1l,m1u)
hypS1<-c(2,200)
# prior on the mode and sample size  of the beta law for the second period
hypMu2<-c(m2l,m2u)
hypS2<-c(2,200)

hyper<-list(mu1=hypMu1,mu2=hypMu2,s1=hypS1,s2=hypS2)


ptm <- proc.time()
# ------------------------
## MCMC settings
# initialization of the algorithm (to be changed if no uniform priors are used)

prodLikInit = 0

while(prodLikInit == 0 | is.na(prodLikInit) | is.infinite(prodLikInit))
{
  meanAdapt<-c(runif(1,hypMu1[1],hypMu1[2]),
               runif(1,hypMu2[1],hypMu2[2]),
               runif(1,hypS1[1],hypS1[2]),
               runif(1,hypS2[1],hypS2[2]))
  
  varAdapt<-(2.38/sqrt(4))*rep(1,4)
  
  lambda<-rep(1,length(meanAdapt))
  state_algo<-new("stateAlgo",probaAcc=numeric(0),
                  lambda=lambda,
                  adaptMean=meanAdapt,
                  adaptVar=varAdapt,
                  accRate=numeric(0))
  
  chain<-data.frame(mu1=numeric(0),mu2=numeric(0),s1=numeric(0),s2=numeric(0))
  ar<-data.frame(ar1=numeric(0),ar2=numeric(0),ar3=numeric(0),ar4=numeric(0),
                 l1=numeric(0),l2=numeric(0),l3=numeric(0),l4=numeric(0),
                 p1=numeric(0),p2=numeric(0),p3=numeric(0),p4=numeric(0))
  ar[1,]<-c(rep(0,4),lambda,rep(0,4))
  chain[1,]<-meanAdapt
  
  # making sure the starting value corresponds to a finite likelihood value
  size_data = nrow(data)
  
  # define a and b param for beta distriutions according to mode and sample size
  # for the candidate
  a1<-chain[1,"mu1"]*(chain[1,"s1"]-2)+1
  b1<-(chain[1,"s1"]-2)*(1-chain[1,"mu1"])+1
  a2<-chain[1,"mu2"]*(chain[1,"s2"]-2)+1
  b2<-(chain[1,"s2"]-2)*(1-chain[1,"mu2"])+1
  
  # we define the two copulas, one for the candidate param and one for the current param
  mat<-as.matrix(data)
  likInit <- dbeta(mat[,1],a1,b1)*dbeta(mat[,2],a1,b1)*dbeta(mat[,3],a2,b2)*dbeta(mat[,4],a2,b2)
  prodLikInit <- prod(likInit)
  print(prodLikInit)
}

# Chain settings
# Type of algorithm : either "hgs" or "mh"
algo<-"hgs" 
# Type of proposal : either "Gl" or "CW"
prop<-"CW" 
# Type of distribution for the proposal : either "unif", "norm" or "beta"
dist<-"norm" 

# run the MCMC algorithm
pb<-txtProgressBar(min=1,max=nb_burn+nb_mcmc,style = 3)
plotInterv<-500
mixProp<-TRUE

tolConv<-0.005
conv<-FALSE
m<-2

# Vector of deviances
dev<-rep(0,nb_mcmc)

var<-rep(NA,4)
bsup<-rep(NA,4)

system.time({
  while(m < (nb_burn+nb_mcmc))
  {
    if(mixProp){
      algo <- "hgs"
      prop <- sample(c("CW","Gl"), size = 1)
      dist <- sample(c("unif","norm"), size = 1)
    }
    
    result<-generateMCMC_indep(chain[m-1,1:4],state_algo,algo,prop,dist,hyper,data)
    chain<-rbind(chain,result$chain)
    state_algo<-adaptAlgo(state_algo,m,unlist(chain[m,]),result$lambda,result$accRate,algo,prop,stochStep = 0.7)
    
    a<-updateAccRate(as.matrix(chain),m,ar[(m-1),])
    l<-state_algo@lambda
    pr<-state_algo@probaAcc
    ar<-rbind(ar,unlist(c(a,l,pr)))
    
    if (m > nb_burn)
    {
      a1<-chain$mu1[m]*(chain$s1[m]-2)+1
      a2<-chain$mu2[m]*(chain$s2[m]-2)+1
      b1<-(chain$s1[m]-2)*(1-chain$mu1[m])+1
      b2<-(chain$s2[m]-2)*(1-chain$mu2[m])+1
      
      mat<-as.matrix(data)
      num<-dbeta(mat[,1],a1,b1)*dbeta(mat[,2],a1,b1)*dbeta(mat[,3],a2,b2)*dbeta(mat[,4],a2,b2)
      dev[m-nb_burn]<-(-2)*sum(log(na.omit(num)))   
      saveRDS(dev,file=paste(landuseCat,"dev_indep",nbrep,sep="_"))       
    }
    
    rownames(chain)<-NULL
    
    saveRDS(chain,file=paste(landuseCat,"chain_indep",nbrep,sep="_"))
    
    m<-m+1
    setTxtProgressBar(pb,m)
  }
})

# reset rownames to avoid very long ones
chain$iter<-seq(1,nrow(chain))
saveRDS(chain,file=paste(landuseCat,"chain_indep",nbrep,sep="_"))

ptmEnd<-proc.time() - ptm
print(ptmEnd)




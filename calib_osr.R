# Floral coverage model

options(echo=FALSE)
args<-commandArgs(trailingOnly = TRUE)
print(args)

args <- c("osr","frank","1","0.5","1","100000","1000")

landuseCat <-args[1]
copFamily <- args[2]
nbrep <- as.integer(args[3])
m1l <- as.numeric(args[4])
m1u <- as.numeric(args[5])
nb_mcmc <- as.numeric(args[6])
nb_burn <- as.numeric(args[7])
rm(args)

setwd("/home/baeyc/Code/es2cropBB/pollination/calibration/floralCoverage/")

# load libraries
library(R.oo)
library(gridExtra)
library(reshape)

library(EnvStats)

source("helpers.R")
source("multiplots.R")

# load data -> we have one matrix of observations per landscape category
source("loadObsData.R")
data<-eval(parse(text=landuseCat))
data[data==1] <- 0.99

## Priors
# 1. Hyperparameters (to be changed for each landscape category)
# prior on the mode and sample size of the beta law for the first period
hypMu1<-c(m1l,m1u)
hypS1<-c(2,100)
# prior on tau parameters used to defined the copula
hypTauW<-c(0,1)

hyper<-list(mu1=hypMu1,s1=hypS1,tauW=hypTauW)

# 2. Choice of the copula family
# one of : gumbel, clayton, frank, joe or amh
if (copFamily == "gumbel"){
  typeHAC<-1
}else if (copFamily == "clayton"){
  typeHAC<-3
}else if (copFamily == "frank"){
  typeHAC<-5
}else if (copFamily == "amh"){
  typeHAC<-9
}else if (copFamily == "joe"){
  typeHAC<-7
}  

ptm <- proc.time()
# ------------------------
## MCMC settings
# initialization of the algorithm (to be changed if no uniform priors are used)

prodLikInit = 0

testedValues <- numeric()
trials <- 0

while((prodLikInit == 0 | is.na(prodLikInit)) & trials < 10000)
{
  meanAdapt<-c(runif(1,hypMu1[1],hypMu1[2]),
               runif(1,hypS1[1],hypS1[2]),
               runif(1,hypTauW[1],hypTauW[2]))
  
  varAdapt<-(2.38/sqrt(3))*rep(1,3)
  
  lambda<-rep(1,length(meanAdapt))
  state_algo<-new("stateAlgo",probaAcc=numeric(0),
                  lambda=lambda,
                  adaptMean=meanAdapt,
                  adaptVar=varAdapt,
                  accRate=numeric(0))
  
  chain<-data.frame(mu1=numeric(0),s1=numeric(0),tauW=numeric(0))
  ar<-data.frame(ar1=numeric(0),ar2=numeric(0),ar3=numeric(0),
                 l1=numeric(0),l2=numeric(0),l3=numeric(0),
                 p1=numeric(0),p2=numeric(0),p3=numeric(0))
  ar[1,]<-c(rep(0,3),lambda,rep(0,3))
  chain[1,]<-meanAdapt
  
  testedValues <- rbind(testedValues,meanAdapt)
  
  # making sure the starting value corresponds to a finite likelihood value
  size_data = nrow(data)
  # we compute theta, the parameters of the copula, according to Kendall's tau values
  if (copFamily == "gumbel"){
    thetaW<-copGumbel@iTau(chain[1,"tauW"])
  }else if (copFamily == "clayton"){
    thetaW<-copClayton@iTau(chain[1,"tauW"])
  }else if (copFamily == "frank"){
    thetaW<-copFrank@iTau(chain[1,"tauW"])
  }else if (copFamily == "joe"){
    thetaW<-copJoe@iTau(chain[1,"tauW"])
  }
  
  # define a and b param for beta distriutions according to mode and sample size
  # for the candidate
  a1<-chain[1,"mu1"]*(chain[1,"s1"]-2)+1
  b1<-(chain[1,"s1"]-2)*(1-chain[1,"mu1"])+1
  #a1 <- chain[1,"mu1"]
  #b1 <- chain[1,"s1"]
  
  # we define the two copulas, one for the candidate param and one for the current param
  mat<-as.matrix(data)
  cop <- hac(type=typeHAC,tree=list("y1","y2",thetaW))
  marginMat <- mat * 0
  marginMat[,1] <- pbeta(mat[,1],a1,b1)
  marginMat[,2] <- pbeta(mat[,2],a1,b1)
  likInit <- dHAC(marginMat,cop)*dlnormTrunc(mat[,1],a1,b1)*dlnormTrunc(mat[,2],a1,b1)
  prodLikInit <- prod(likInit)
  
  trials <- trials+1
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

var<-rep(NA,3)
bsup<-rep(NA,3)

system.time({
  while(m < (nb_burn+nb_mcmc))
  {
    if(mixProp){
      algo <- "hgs"
      prop <- sample(c("CW","Gl"), size = 1)
      dist <- sample(c("unif","norm"), size = 1)
    }
    
    result<-generateMCMC_OSR(chain[m-1,1:3],state_algo,algo,prop,dist,hyper,copFamily,data)
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
      
      if (copFamily == "gumbel"){
        theta1<-copGumbel@iTau(chain$tauW[m])
        typeHAC<-1
      }else if (copFamily == "clayton"){
        theta1<-copClayton@iTau(chain$tauW[m])
        typeHAC<-3
      }else if (copFamily == "frank"){
        theta1<-copFrank@iTau(chain$tauW[m])
        typeHAC<-5
      }else if (copFamily == "amh"){
        theta1<-copAMH@iTau(chain$tauW[m])
        typeHAC<-9
      }else if (copFamily == "joe"){
        theta1<-copJoe@iTau(chain$tauW[m])
        typeHAC<-7
      }      
      
      mat<-as.matrix(data)
      cop<-hac(type=typeHAC,tree=list("y1","y2",theta1))
      num<-dHAC(mat,cop)*dbeta(mat[,1],a1,b1)*dbeta(mat[,2],a1,b1)
      dev[m-nb_burn]<-(-2)*sum(log(na.omit(num)))   
      saveRDS(dev,file=paste(landuseCat,"dev",copFamily,nbrep,sep="_"))       
    }
    
    rownames(chain)<-NULL
    
    saveRDS(chain,file=paste(landuseCat,"chain",copFamily,nbrep,sep="_"))
    
    m<-m+1
    setTxtProgressBar(pb,m)
  }
})

# reset rownames to avoid very long ones
chain$iter<-seq(1,nrow(chain))
saveRDS(chain,file=paste(landuseCat,"chain",copFamily,nbrep,sep="_"))

ptmEnd<-proc.time() - ptm
print(ptmEnd)




# Floral coverage model

#setwd("/lunarc/nobackup/users/baeyc/floralCoverage/")

lu <- "pasture"
nb_mcmc <- 100
nb_burn <- 10
nbrep<-1 # nb of replicates of the mcmc algo
# Type of algorithm : either "hgs" or "mh"
algo<-"mh" 
# Type of proposal : either "Gl" or "CW"
prop<-"Gl" 
# Type of distribution for the proposal : either "unif" or "norm"
dist<-"unif" 
mixProp<-FALSE
eval(parse(text=paste(commandArgs(trailingOnly = TRUE), collapse=";")))

# load libraries
library(R.oo)
library(ggplot2)
library(gridExtra)
library(reshape)

source("helpers.R")
source("multiplots.R")

# load data -> we have one matrix of observations per landscape category
source("loadObsData.R")

## First category -> pasture
if (lu == "pasture")
{
  hypMu1<-c(0,0.2)
  hypMu2<-c(0,0.2)
} else if (lu == "fieldEdge")
{
  hypMu1<-c(0,0.2)
  hypMu2<-c(0,0.2)
} else if (lu == "flowStrip")
{
  hypMu1<-c(0,0.3)
  hypMu2<-c(0,0.3)
}

## Priors
# 1. Hyperparameters (to be changed for each landscape category)
# prior on the mode and sample size of the beta law for the first period
# prior on the mode and sample size  of the beta law for the second period
hypS1<-c(2,200)
hypS2<-c(2,200)
# prior on tau parameters used to defined the copula
hypTauW<-c(0,1)
hypTauB<-0 # only the lower bound, as the upper bound will be given by hypTauW

hyper<-list(mu1=hypMu1,mu2=hypMu2,s1=hypS1,s2=hypS2,tauW=hypTauW,tauB=hypTauB)

# 2. Choice of the copula family
# one of : gumbel, clayton, frank, joe or amh
copFamily<-"joe"
typeHAC<-3

for (rep in 1:nbrep){
  ptm <- proc.time()
  # ------------------------
  ## MCMC settings
  # initialization of the algorithm (to be changed if no uniform priors are used)
  #meanAdapt<-c(mean(hypMu1),mean(hypMu2),mean(hypS1),mean(hypS2),mean(hypTauW),mean(c(hypTauB,hypTauW)))
  meanAdapt<-c(runif(1,hypMu1[1],hypMu1[2]),
               runif(1,hypMu2[1],hypMu2[2]),
               runif(1,hypS1[1],hypS1[2]),
               runif(1,hypS2[1],hypS2[2]),
               runif(1,hypTauW[1],hypTauW[2]))
  meanAdapt<-c(meanAdapt,runif(1,0,meanAdapt[5]))
  
  #varAdapt<-c((hypMu1[1]-hypMu1[2])^2/12,
  #           (hypMu2[1]-hypMu2[2])^2/12,
  #          (hypS1[1]-hypS1[2])^2/12,
  #         (hypS2[1]-hypS2[2])^2/12,
  #        (hypTauW[1]-hypTauW[2])^2/12,
  #       (hypTauB-mean(hypTauW))^2/12)
  varAdapt<-(2.38/sqrt(6))*rep(1,6)
  #varAdapt<-meanAdapt
  lambda<-rep(1,length(meanAdapt))
  state_algo<-new("stateAlgo",probaAcc=numeric(0),
                  lambda=lambda,
                  adaptMean=meanAdapt,
                  adaptVar=varAdapt,
                  accRate=numeric(0))
  
  chain<-data.frame(mu1=numeric(0),mu2=numeric(0),s1=numeric(0),s2=numeric(0),tauW=numeric(0),tauB=numeric(0))
  ar<-data.frame(ar1=numeric(0),ar2=numeric(0),ar3=numeric(0),ar4=numeric(0),ar5=numeric(0),ar6=numeric(0),
                 l1=numeric(0),l2=numeric(0),l3=numeric(0),l4=numeric(0),l5=numeric(0),l6=numeric(0),
                 p1=numeric(0),p2=numeric(0),p3=numeric(0),p4=numeric(0),p5=numeric(0),p6=numeric(0))
  ar[1,]<-c(rep(0,6),lambda,rep(0,6))
  chain[1,]<-meanAdapt
  
  # run the MCMC algorithm
  pb<-txtProgressBar(min=1,max=nb_burn+nb_mcmc,style = 3)
  plotInterv<-50
  
  tolConv<-0.005
  conv<-FALSE
  m<-2
  
  # Defining layout for monitoring graphs
  layout(matrix(c(1,1,1,1,1,1,2,3,4,5,6,7),nr=2))
  
  # Vector of deviances
  dev<-rep(0,nb_mcmc)
  
  var<-rep(NA,6)
  bsup<-rep(NA,6)
  
  try(
    while((m < nb_burn) || !conv)
    {
      if(mixProp){
        algo<-sample(c("mh","hgs"),size = 1)
        prop<-sample(c("CW","Gl"),size = 1)
      }
      
      result<-generateMCMC(chain[m-1,1:6],state_algo,algo,prop,dist,hyper,copFamily,pasture)
      chain<-rbind(chain,result$chain)
      state_algo<-adaptAlgo(state_algo,m,unlist(chain[m,]),result$lambda,result$accRate,algo,prop,0.7)
      
      a<-updateAccRate(as.matrix(chain),m,ar[(m-1),])
      l<-state_algo@lambda
      pr<-state_algo@probaAcc
      print(state_algo@adaptVar)
      ar<-rbind(ar,unlist(c(a,l,pr)))
      
      if (m > nb_burn)
      {
        a1<-chain$mu1[m]*(chain$s1[m]-2)+1
        a2<-chain$mu2[m]*(chain$s2[m]-2)+1
        b1<-(chain$s1[m]-2)*(1-chain$mu1[m])+1
        b2<-(chain$s2[m]-2)*(1-chain$mu2[m])+1
        
        if (copFamily == "gumbel"){
          theta1<-copGumbel@iTau(chain$tauW[m])
          theta2<-copGumbel@iTau(chain$tauB[m])
          typeHAC<-1
        }else if (copFamily == "clayton"){
          theta1<-copClayton@iTau(chain$tauW[m])
          theta2<-copClayton@iTau(chain$tauB[m])
          typeHAC<-3
        }else if (copFamily == "frank"){
          theta1<-copFrank@iTau(chain$tauW[m])
          theta2<-copFrank@iTau(chain$tauB[m]) 
          typeHAC<-5
        }else if (copFamily == "amh"){
          theta1<-copAMH@iTau(chain$tauW[m])
          theta2<-copAMH@iTau(chain$tauB[m])
          typeHAC<-9
        }else if (copFamily == "joe"){
          theta1<-copJoe@iTau(chain$tauW[m])
          theta2<-copJoe@iTau(chain$tauB[m])   
          typeHAC<-7
        }      
        
        mat<-as.matrix(pasture)
        cop<-hac(type=typeHAC,tree=list(list("y1","y2",theta1),list("y3","y4",theta1),theta2))
        num<-dHAC(mat,cop)*dbeta(mat[,1],a1,b1)*dbeta(mat[,2],a1,b1)*dbeta(mat[,3],a2,b2)*dbeta(mat[,4],a2,b2)
        dev[m-nb_burn]<-(-2)*sum(log(na.omit(num)))
        
        if ((m-nb_burn) %% plotInterv == 0)
        {
          var<-rbind(var,sapply(chain[(nb_burn+1):m,1:6],FUN=batchVar))
          k<-nrow(var)
          bsup<-rbind(bsup,qt(1-0.025,df=1)*sqrt(unlist(var[k,])/(m-nb_burn)))    
          #binf<-rbind(binf,qt(0.025,df=1)*sqrt(unlist(var[k,])/(m-nb_burn)))  
          
          #par(new=FALSE)
          #for(i in 1:6){
          #  plot(seq(1:k),var[,i],type="l",ylim=c(min(na.omit(var[,i]+bsup[,i])*0.9),max(na.omit(var[,i]+bsup[,i])*1.1)),xlim=c(0,max(m,nb_burn+nb_mcmc)/plotInterv),xlab="",ylab="")
          #  par(new=TRUE)
          #  plot(seq(1:k),var[,i]+bsup[,i],type="l",lty=2,ylim=c(min(na.omit(var[,i]+bsup[,i])*0.9),max(na.omit(var[,i]+bsup[,i])*1.1)),xlim=c(0,max(m,nb_burn+nb_mcmc)/plotInterv),xlab="",ylab="")
          #}
          
          if (max(bsup[k,]/var[k,]) < tolConv || m == (nb_mcmc+nb_burn))
            conv<-TRUE
        }      
      }
      
      m<-m+1
      setTxtProgressBar(pb,m)
    }
  )
  
  # reset rownames to avoid very long ones
  rownames(ar)<-NULL
  rownames(chain)<-NULL
  rownames(bsup)<-NULL
  rownames(var)<-NULL
  
  chain$iter<-seq(1,nrow(chain))
  ar$iter<-seq(1,nrow(ar))
  batch<-data.frame(cbind(var,var+bsup))
  batch$iter<-seq(1,nrow(batch))
  
  #Save results
  saveRDS(chain,file=paste("chain",copFamily,rep,".rds",sep="_"))
  saveRDS(meanAdapt,file=paste("init",copFamily,rep,".rds",sep="_"))
  saveRDS(ar,file=paste("AR",copFamily,rep,".rds",sep="_"))
  saveRDS(batch,file=paste("batch",copFamily,rep,".rds",sep="_"))
  saveRDS(dev,file=paste("dev",copFamily,rep,".rds",sep="_"))
  
  ptmEnd<-proc.time() - ptm
  print(ptmEnd)
  
}# end rep chain

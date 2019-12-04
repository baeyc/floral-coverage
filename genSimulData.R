# Load libraries
library(copula)
library(HAC)



# Parameters
# 1. Hyperparameters (to be changed for each landscape category)
# prior on the mode and sample size of the beta law for the first period
hypMu1<-c(0.05,0.10)
hypS1<-c(20,30)
# prior on the mode and sample size  of the beta law for the second period
hypMu2<-c(0.05,0.10)
hypS2<-c(20,30)
# prior on tau parameters used to defined the copula
hypTauW<-c(0.6,0.8)
hypTauB<-c(0.4,0.6)

# Generate parameters from priors distributions
mu1<-runif(1,min=hypMu1[1],max=hypMu1[2])
mu2<-runif(1,min=hypMu2[1],max=hypMu2[2])

tauB<-runif(1,hypTauB[1],hypTauB[2])
tauW<-runif(1,hypTauW[1],hypTauW[2])

thetaW<-copGumbel@iTau(tauW)
thetaB<-copGumbel@iTau(tauB)

# we define the two copulas, one for the candidate param and one for the current param
copCand<-hac(type=1,tree=list(list("y1","y2",thetaW),list("y3","y4",thetaW),thetaW))

rcop <- rHAC(50,copCand)
plot(rcop)

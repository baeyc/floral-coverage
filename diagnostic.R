library(ggplot2)
library(gridExtra)
library(reshape)
library(HAC)
library(coda)

setwd("/home/baeyc/Code/es2cropBB/pollination/calibration/floralCoverage/")

source("multiplots.R")

# data
source("loadObsData.R")

## First category -> pasture
pasture<-pastMat[,c(3,4,5,6)]
pasture<-pasture/100
names(pasture)<-c("y1","y2","y3","y4")

# Best model and family -> Gumbel 2
chainPast <- readRDS("Results/chain_gumbel_2")
chainPast$ar <- rep(0,nrow(chainPast))
chainPast$ar[1] <- 0
for (i in 2:nrow(chainPast))
{
  if (chainPast$mu1[i] != chainPast$mu1[i-1])
    chainPast$ar[i] = chainPast$ar[i-1]+1
  else
    chainPast$ar[i] = chainPast$ar[i-1]
}
chainPast$iter <- seq(1,nrow(chainPast),1)
chainPast$ar2 <- chainPast$ar/chainPast$iter

par(mfrow=c(1,1))
plot(chainPast$iter,chainPast$ar2,type="l",xlim=c(0,100000))
abline(h=0.234,lty=2)

codaPast <- mcmc(chainPast, start=1000)

# Acceptance rate
1-rejectionRate(codaPast)

# Autocorrelation
autocorr.plot(codaPast)
autocorr.diag(codaPast)

# Cumulative quantile plots
thinnedPast <- window(codaPast, thin = 100)
cumuplot(thinnedPast, probs=c(0.025,0.5,0.975))

# Diagnostics
heidel.diag(codaPast)

chainPast2 <- readRDS("Results/chain_gumbel_3")
codaPast2 <- mcmc(chainPast2, start = 1000)
listPast <- mcmc.list(codaPast,codaPast2)
gelman.diag(listPast)
gelman.plot(listPast)

## Field edges
fieldEdge<-feMat[,c(3,4,5,6)]
fieldEdge<-fieldEdge/100
names(fieldEdge)<-c("y1","y2","y3","y4")
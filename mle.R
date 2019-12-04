# Floral coverage model
# Frequentist estimation

setwd("/lunarc/nobackup/users/baeyc/floralCoverage/")

# load libraries
library(R.oo)
library(ggplot2)
library(gridExtra)
library(reshape)
library(HAC)

source("multiplots.R")

# load data -> we have one matrix of observations per landscape category
source("loadObsData.R")

## First category -> pasture
pasture<-pastMat[,c(3,4,5,6)]
pasture<-pasture/100
names(pasture)<-c("y1","y2","y3","y4")

copCand<-hac(type=1,tree=list(list("y1","y2",0.8),list("y3","y4",0.8),0.6))
plot(copCand)

estCOP<-estimate.copula(na.omit(pasture),margins = "edf")
plot(estCOP)

tau1<-copGumbel@tau(7.38)
tau2<-copGumbel@tau(5.07)
tau3<-copGumbel@tau(4.39)

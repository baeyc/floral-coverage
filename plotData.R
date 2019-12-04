# Floral coverage model

library(copula)
library(GGally)
library(ggplot2)
library(extrafont)
library(fontcm)

Sys.setenv(R_GSCMD = "C:/Program Files (x86)/gs/gs9.16/bin/gswin32c.exe")

# simulations
theta0 <- copClayton@iTau(0.5)
theta1 <- copClayton@iTau(0.8)

copC<-onacopula("Clayton",C(theta0, NULL,list(C(theta1,c(1,2)),C(theta1,c(3,4)))))
#copG<-onacopula("Gumbel",C(theta0, NULL,list(C(theta1,c(1,2)),C(theta1,c(3,4)))))

n<-16
UC<-rnacopula(n,copC)
#UG<-rnacopula(n,copG)

mode<-0.45
sampleSize<-50
a<-mode*(sampleSize-2)+1
b<-(sampleSize-2)*(1-mode)+1
YC<-qbeta(UC,shape1 = a,shape2 = b)
#YG<-qbeta(UG,shape1 = a,shape2 = b)

ggpairs(YC,columnLabels = c("Y1","Y2","Y3","Y4"))

# data
setwd("C:/Users/Charlotte/Code/es2crop/empirical/STEP/")

data<-read.table("STEPsitesWithLUCat.txt",header=TRUE,sep="\t")
data$id<-paste(data$Lokal,data$Year,sep="_")
data<-data[,c(5,11,12)]

pasture<-data[data$landCategory==1,]
fieldEdge<-data[data$landCategory==2,]
osr<-data[data$landCategory==3,]
meadow<-data[data$landCategory==4,]
flowStrip<-data[data$landCategory==5,]
fieldEdgeFl<-data[data$landCategory==6,]

pasture$time<-c(1,2,3,4)
fieldEdge$time<-c(1,2,3,4)
meadow$time<-c(1,2,3,4)
flowStrip$time<-c(1,2,3,4)
fieldEdgeFl$time<-c(1,2,3,4)
osr$time<-c(1,2)

# Plotting data -> reshape dataframes into matrices
pastMat<-reshape(pasture,v.names="TransectFlowerCover",timevar = "time",idvar = "id",direction = "wide")
feMat<-reshape(fieldEdge,v.names="TransectFlowerCover",timevar = "time",idvar = "id",direction = "wide")
osrMat<-reshape(osr,v.names="TransectFlowerCover",timevar = "time",idvar = "id",direction = "wide")
meadMat<-reshape(meadow,v.names="TransectFlowerCover",timevar = "time",idvar = "id",direction = "wide")
fSMat<-reshape(flowStrip,v.names="TransectFlowerCover",timevar = "time",idvar = "id",direction = "wide")
feFMat<-reshape(fieldEdgeFl,v.names="TransectFlowerCover",timevar = "time",idvar = "id",direction = "wide")

pdf("PastureData.pdf",fonts="Garamond")
ggpairs(pastMat[,c(3,4,5,6)],columnLabels = c("Date 1","Date 2","Date 3","Date 4"))
dev.off()
embed_fonts("PastureData.pdf")

pdf("FieldEdgeData.pdf",fonts="Garamond")
ggpairs(feMat[,c(3,4,5,6)],columnLabels = c("Date 1","Date 2","Date 3","Date 4"))
dev.off()
embed_fonts("FieldEdgeData.pdf")

pdf("OSRData.pdf",fonts="Garamond")
ggpairs(osrMat[,c(3,4)],columnLabels = c("Date 1","Date 2"))
dev.off()
embed_fonts("OSRData.pdf")

pdf("MeadowData.pdf",fonts="Garamond")
ggpairs(meadMat[,c(3,4,5,6)],columnLabels = c("Date 1","Date 2","Date 3","Date 4"))
dev.off()
embed_fonts("MeadowData.pdf")

pdf("FlowStripsData.pdf",fonts="Garamond")
ggpairs(fSMat[,c(3,4,5,6)],columnLabels = c("Date 1","Date 2","Date 3","Date 4"))
dev.off()
embed_fonts("FlowStripsData.pdf")

pdf("FieldEdgeFlData.pdf",fonts="Garamond")
ggpairs(feFMat[,c(3,4,5,6)],columnLabels = c("Date 1","Date 2","Date 3","Date 4"))
dev.off()
embed_fonts("FieldEdgeFlData.pdf")


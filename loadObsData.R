# Load real data of floral coverage

curWD<-getwd()

data<-read.table("STEPsitesWithLUCat.txt",header=TRUE,sep="\t")
data$id<-paste(data$Lokal,data$Year,sep="_")
data<-data[,c(5,11,12)]

pasture<-data[data$landCategory %in% c(1,4),]
fieldEdge<-data[data$landCategory==2,]
osr<-data[data$landCategory==3,]
#meadow<-data[data$landCategory==4,]
flowStrip<-data[data$landCategory==5,]
fieldEdgeFl<-data[data$landCategory==6,]

pasture$time<-c(1,2,3,4)
fieldEdge$time<-c(1,2,3,4)
#meadow$time<-c(1,2,3,4)
flowStrip$time<-c(1,2,3,4)
fieldEdgeFl$time<-c(1,2,3,4)
osr$time<-c(1,2)

pastMat<-reshape(pasture,v.names="TransectFlowerCover",timevar = "time",idvar = "id",direction = "wide")
feMat<-reshape(fieldEdge,v.names="TransectFlowerCover",timevar = "time",idvar = "id",direction = "wide")
osrMat<-reshape(osr,v.names="TransectFlowerCover",timevar = "time",idvar = "id",direction = "wide")
#meadMat<-reshape(meadow,v.names="TransectFlowerCover",timevar = "time",idvar = "id",direction = "wide")
fSMat<-reshape(flowStrip,v.names="TransectFlowerCover",timevar = "time",idvar = "id",direction = "wide")
feFMat<-reshape(fieldEdgeFl,v.names="TransectFlowerCover",timevar = "time",idvar = "id",direction = "wide")

setwd(curWD)


## Give column names to data matrix to be used in copula functions
pasture<-pastMat[,c(3,4,5,6)]
pasture<-pasture/100
names(pasture)<-c("y1","y2","y3","y4")

fieldEdge<-feMat[,c(3,4,5,6)]
fieldEdge<-fieldEdge/100
names(fieldEdge)<-c("y1","y2","y3","y4")

osr<-osrMat[,c(3,4)]
osr<-osr/100
names(osr)<-c("y1","y2")

#meadow<-meadMat[,c(3,4,5,6)]
#meadow<-meadow/100
#names(meadow)<-c("y1","y2","y3","y4")

flowStrip<-fSMat[,c(3,4,5,6)]
flowStrip<-flowStrip/100
names(flowStrip)<-c("y1","y2","y3","y4")

fieldEdgeFl<-feFMat[,c(3,4,5,6)]
fieldEdgeFl<-fieldEdgeFl/100
names(fieldEdgeFl)<-c("y1","y2","y3","y4")



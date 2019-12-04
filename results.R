
setwd("/home/baeyc/Code/es2cropBB/pollination/calibration/floralCoverage/Results/")
library(copula)
library(HAC)

source("../loadObsData.R")

## -------------------------------------------------
##                     Gumbel
## -------------------------------------------------

init1<-readRDS("InitGumbel1")
init2<-readRDS("IniGumbelt2")
init3<-readRDS("InitGumbel3")

batch1<-readRDS("BatchGumbel1")
batch2<-readRDS("BatchGumbel2")
batch3<-readRDS("BatchGumbel3")

chain1<-readRDS("chain_gumbel_1")
chain2<-readRDS("chain_gumbel_2")
chain3<-readRDS("chain_gumbel_3")
chain1$iter<-seq(1,nrow(chain1))
chain2$iter<-seq(1,nrow(chain2))
chain3$iter<-seq(1,nrow(chain3))
#write.table(chain1,"Chain1.txt",row.names = FALSE)
#write.table(chain2,"Chain2.txt",row.names = FALSE)
#write.table(chain3,"Chain3.txt",row.names = FALSE)

chains<-cbind(c(chain1$mu1[(nb_burn+1):100000],rep(NA,100000)),chain2$mu1[(nb_burn+1):200000],chain3$mu1[(nb_burn+1):200000])
chainsMu1<-melt(chains,id="V4")
chains<-cbind(c(chain1$mu2[(nb_burn+1):100000],rep(NA,100000)),chain2$mu2[(nb_burn+1):200000],chain3$mu2[(nb_burn+1):200000])
chainsMu2<-melt(chains,id="V4")
chains<-cbind(c(chain1$s1[(nb_burn+1):100000],rep(NA,100000)),chain2$s1[(nb_burn+1):200000],chain3$s1[(nb_burn+1):200000])
chainss1<-melt(chains,id="V4")
chains<-cbind(c(chain1$s2[(nb_burn+1):100000],rep(NA,100000)),chain2$s2[(nb_burn+1):200000],chain3$s2[(nb_burn+1):200000])
chainss2<-melt(chains,id="V4")
chains<-cbind(c(chain1$tauW[(nb_burn+1):100000],rep(NA,100000)),chain2$tauW[(nb_burn+1):200000],chain3$tauW[(nb_burn+1):200000])
chainstauW<-melt(chains,id="V4")
chains<-cbind(c(chain1$tauB[(nb_burn+1):100000],rep(NA,100000)),chain2$tauB[(nb_burn+1):200000],chain3$tauB[(nb_burn+1):200000])
chainstauB<-melt(chains,id="V4")

# Criteria computation
# Gambel copula (type 1 in hac package)
dev1<-readRDS("dev_gumbel_1")
dev2<-readRDS("dev_gumbel_2")
dev3<-readRDS("dev_gumbel_3")

# deviance for posterior mean
nb_burn<-1000
mat<-as.matrix(pasture)

a1<-mean(chain1$mu1[(nb_burn+1):nrow(chain1)])*(mean(chain1$s1[(nb_burn+1):nrow(chain1)])-2)+1
a2<-mean(chain1$mu2[(nb_burn+1):nrow(chain1)])*(mean(chain1$s2[(nb_burn+1):nrow(chain1)])-2)+1
b1<-(mean(chain1$s1[(nb_burn+1):nrow(chain1)])-2)*(1-mean(chain1$mu1[(nb_burn+1):nrow(chain1)]))+1
b2<-(mean(chain1$s2[(nb_burn+1):nrow(chain1)])-2)*(1-mean(chain1$mu2[(nb_burn+1):nrow(chain1)]))+1
cop<-hac(type=1,tree=list(list("y1","y2",mean(chain1$tauW[(nb_burn+1):nrow(chain1)])),list("y3","y4",mean(chain1$tauW[(nb_burn+1):nrow(chain1)])),mean(chain1$tauB[(nb_burn+1):nrow(chain1)])))
num<-dHAC(mat,cop)*dbeta(mat[,1],a1,b1)*dbeta(mat[,2],a1,b1)*dbeta(mat[,3],a2,b2)*dbeta(mat[,4],a2,b2)
devHat1<-(-2)*sum(log(na.omit(num)))
###
a1<-mean(chain2$mu1[(nb_burn+1):nrow(chain2)])*(mean(chain2$s1[(nb_burn+1):nrow(chain2)])-2)+1
a2<-mean(chain2$mu2[(nb_burn+1):nrow(chain2)])*(mean(chain2$s2[(nb_burn+1):nrow(chain2)])-2)+1
b1<-(mean(chain2$s1[(nb_burn+1):nrow(chain2)])-2)*(1-mean(chain2$mu1[(nb_burn+1):nrow(chain2)]))+1
b2<-(mean(chain2$s2[(nb_burn+1):nrow(chain2)])-2)*(1-mean(chain2$mu2[(nb_burn+1):nrow(chain2)]))+1
cop<-hac(type=1,tree=list(list("y1","y2",mean(chain2$tauW[(nb_burn+1):nrow(chain2)])),list("y3","y4",mean(chain2$tauW[(nb_burn+1):nrow(chain2)])),mean(chain2$tauB[(nb_burn+1):nrow(chain2)])))
num<-dHAC(mat,cop)*dbeta(mat[,1],a1,b1)*dbeta(mat[,2],a1,b1)*dbeta(mat[,3],a2,b2)*dbeta(mat[,4],a2,b2)
devHat2<-(-2)*sum(log(na.omit(num)))
###
a1<-mean(chain3$mu1[(nb_burn+1):nrow(chain3)])*(mean(chain3$s1[(nb_burn+1):nrow(chain3)])-2)+1
a2<-mean(chain3$mu2[(nb_burn+1):nrow(chain3)])*(mean(chain3$s2[(nb_burn+1):nrow(chain3)])-2)+1
b1<-(mean(chain3$s1[(nb_burn+1):nrow(chain3)])-2)*(1-mean(chain3$mu1[(nb_burn+1):nrow(chain3)]))+1
b2<-(mean(chain3$s2[(nb_burn+1):nrow(chain3)])-2)*(1-mean(chain3$mu2[(nb_burn+1):nrow(chain3)]))+1
cop<-hac(type=1,tree=list(list("y1","y2",mean(chain3$tauW[(nb_burn+1):nrow(chain3)])),list("y3","y4",mean(chain3$tauW[(nb_burn+1):nrow(chain3)])),mean(chain3$tauB[(nb_burn+1):nrow(chain3)])))
num<-dHAC(mat,cop)*dbeta(mat[,1],a1,b1)*dbeta(mat[,2],a1,b1)*dbeta(mat[,3],a2,b2)*dbeta(mat[,4],a2,b2)
devHat3<-(-2)*sum(log(na.omit(num)))

pD<-mean(dev1)-devHat1
DIC1<-devHat1+2*pD # -2218.83015
pD<-mean(dev2)-devHat2
DIC2<-devHat2+2*pD # -2217.24761
pD<-mean(dev3)-devHat3
DIC3<-devHat3+2*pD # -2216.40938


# Histograms posteriors
library(easyGgplot2)
h1<-ggplot2.histogram(data=chainsMu1, xName='value',groupName='X2', backgroundColor="white",removePanelGrid=TRUE,removePanelBorder=TRUE,legendPosition="top",alpha=0.5, position="dodge", scale = "density", show_guide=FALSE)+labs(title=expression(mu[1]),x="",y="")+ theme(plot.title = element_text(size = rel(2)))
h2<-ggplot2.histogram(data=chainsMu2, xName='value',groupName='X2',addMeanLine = FALSE, backgroundColor="white",removePanelGrid=TRUE,removePanelBorder=TRUE, legendPosition="top",alpha=0.5, position="dodge", scale = "density", show_guide=FALSE)+labs(title=expression(mu[2]),x="",y="")+ theme(plot.title = element_text(size = rel(2)))
h3<-ggplot2.histogram(data=chainss1, xName='value',groupName='X2',addMeanLine = FALSE, backgroundColor="white",removePanelGrid=TRUE,removePanelBorder=TRUE, legendPosition="top",alpha=0.5, position="dodge", scale = "density", show_guide=FALSE)+labs(title=expression(s[1]),x="",y="")+ theme(plot.title = element_text(size = rel(2)))
h4<-ggplot2.histogram(data=chainss2, xName='value',groupName='X2',addMeanLine = FALSE, backgroundColor="white",removePanelGrid=TRUE,removePanelBorder=TRUE, legendPosition="top",alpha=0.5, position="dodge", scale = "density", show_guide=FALSE)+labs(title=expression(s[2]),x="",y="")+ theme(plot.title = element_text(size = rel(2)))
h5<-ggplot2.histogram(data=chainstauW, xName='value',groupName='X2',addMeanLine = FALSE, backgroundColor="white",removePanelGrid=TRUE,removePanelBorder=TRUE, legendPosition="top",alpha=0.5, position="dodge", scale = "density", show_guide=FALSE)+labs(title=expression(tau[W]),x="",y="")+ theme(plot.title = element_text(size = rel(2)))
h6<-ggplot2.histogram(data=chainstauB, xName='value',groupName='X2',addMeanLine = FALSE, backgroundColor="white",removePanelGrid=TRUE,removePanelBorder=TRUE, legendPosition="top",alpha=0.5, position="dodge", scale = "density", show_guide=FALSE)+labs(title=expression(tau[B]),x="",y="")+ theme(plot.title = element_text(size = rel(2)))

pdf("Posteriors.pdf",height=10,width=15)
multiplot(h1,h2,h3,h4,h5,h6,cols=3)
dev.off()

# MCMC traces
p1<-ggplot(data=chainsMu1,aes(x=X1,y=value,colour=factor(X2))) + geom_line(size=0.75) + scale_colour_discrete()#name="",values=c("black","darkgrey","darkgrey")) + guides(color="none")
p2<-ggplot(data=chainsMu2,aes(x=X1,y=value,colour=factor(X2))) + geom_line(size=0.75) + scale_colour_discrete()
p3<-ggplot(data=chainss1,aes(x=X1,y=value,colour=factor(X2))) + geom_line(size=0.75) + scale_colour_discrete()
p4<-ggplot(data=chainss2,aes(x=X1,y=value,colour=factor(X2))) + geom_line(size=0.75) + scale_colour_discrete()
p5<-ggplot(data=chainstauW,aes(x=X1,y=value,colour=factor(X2))) + geom_line(size=0.75) + scale_colour_discrete()
p6<-ggplot(data=chainstauB,aes(x=X1,y=value,colour=factor(X2))) + geom_line(size=0.75) + scale_colour_discrete()

pdf("Traces.pdf",height=10,width=15)
multiplot(p1,p2,p3,p4,p5,p6,cols=3)
dev.off()
# 
# ## Graphes priors/posteriors
# # mu1
# mu<-mean(median(chain1$mu1),median(chain2$mu1),median(chain3$mu1))
# s<-mean(median(chain1$s1),median(chain2$s1),median(chain3$s1))
# a1<-mu*(s-2)+1
# b1<-(s-2)*(1-mu)+1
# 
# cairo_pdf("prior.pdf",width=8,height = 7)
# plot(seq(0,0.5,0.001),dunif(seq(0,0.5,0.001),0,0.2),type="l",ylim=c(0,220),xlim=c(0,0.25),col="red",lwd=3,xlab="",ylab="",cex.axis=1.5)
# legend(0.15,220,legend=c("Prior","Posterior","Floral coverage"),lwd=3,lty=1,col=c("red","darkgreen","orange"),bty="n",cex=1.5)
# legend(0.167,175,legend=c("Observations","Observations"),pch="|",col=c("black","grey"),bty="n",cex=1.5)
# dev.off()
# cairo_pdf("priorObs.pdf",width=8,height = 7)
# points(pasture$y1,rep(-5,length(pasture$y1)),pch="|")
# points(pasture$y2,rep(-5,length(pasture$y2)),pch="|",col="darkgrey")
# cairo_pdf("posteriorMode.pdf",width=8,height = 7)
# par(new=TRUE)
# plot(density(chain1$mu1[(nb_burn+1):100000]),xlim=c(0,0.25),ylim=c(0,220),lwd=3,col="darkgreen",xlab="",ylab="",main="",axes=FALSE)
# abline(v=median(chain1$mu1[(nb_burn+1):100000]),lwd=3,col="darkgreen",lty=2)
# cairo_pdf("floralCovMedaPost.pdf",width=8,height = 7)
# par(new=TRUE)
# plot(seq(0,1,0.001),dbeta(seq(0,1,0.001),a1,b1),ylim=c(0,220),xlim=c(0,0.25),type="l",axes=FALSE,xlab="",ylab="",lwd=3,col="orange")
# legend(0.19,220,legend=c("Prior","Posterior","Floral coverage"),lwd=3,lty=1,col=c("red","darkgreen","orange"),bty="n")
# legend(0.201,190,legend=c("Observations","Observations"),pch="|",col=c("black","grey"),bty="n")
# 
# 
# library(triangle)
# 
# cairo_pdf("Pasture_1stPeriod.pdf",height = 8,width=7)
# plot(seq(0,1,0.001),dtriangle(seq(0,1,0.001),a = 0.03,b=0.25,c=0.07),cex.axis=1.5,cex.lab=1.5,lwd=3,type="l",ylim=c(0,35),xlim=c(0,0.4),xlab="% floral coverage",ylab="Density")
# points(pasture$y1,rep(-0.8,length(pasture$y1)),pch="|")
# points(pasture$y2,rep(-0.8,length(pasture$y2)),pch="|",col="grey")
# par(new=TRUE)
# plot(seq(0,1,0.001),dbeta(seq(0,1,0.001),a1,b1),ylim=c(0,35),xlim=c(0,0.4),type="l",axes=FALSE,xlab="",ylab="",lwd=3,col="orange")
# legend(0.25,35,legend=c("Triang. dist","New dist."),col=c("black","orange"),bty="n",cex=1.5,lwd=3)
# legend(0.285,31,legend=c("obs","obs"),col=c("black","grey"),pch="|",bty="n",cex=1.5)
# dev.off()
# 
# mu<-mean(median(chain1$mu2),median(chain2$mu2),median(chain3$mu2))
# s<-mean(median(chain1$s2),median(chain2$s2),median(chain3$s2))
# a2<-mu*(s-2)+1
# b2<-(s-2)*(1-mu)+1
# 
# cairo_pdf("Pasture_2ndPeriod.pdf",height = 8,width=7)
# plot(seq(0,1,0.001),dtriangle(seq(0,1,0.001),a = 0.01,b=0.10,c=0.03),cex.axis=1.5,cex.lab=1.5,lwd=3,type="l",ylim=c(0,70),xlim=c(0,0.4),xlab="% floral coverage",ylab="Density")
# points(pasture$y3,rep(-0.8,length(pasture$y1)),pch="|")
# points(pasture$y4,rep(-0.8,length(pasture$y2)),pch="|",col="grey")
# par(new=TRUE)
# plot(seq(0,1,0.001),dbeta(seq(0,1,0.001),a2,b2),ylim=c(0,70),xlim=c(0,0.4),type="l",axes=FALSE,xlab="",ylab="",lwd=3,col="orange")
# legend(0.25,70,legend=c("Triang. dist","New dist."),col=c("black","orange"),bty="n",cex=1.5,lwd=3)
# legend(0.285,62.5,legend=c("obs","obs"),col=c("black","grey"),pch="|",bty="n",cex=1.5)
# dev.off()
# 
# plot(seq(0,1,0.001),dbeta(seq(0,1,0.001),a1,b1),ylim=c(0,70),xlim=c(0,0.4),type="l")
# par(new=TRUE)
# plot(seq(0,1,0.001),dbeta(seq(0,1,0.001),a2,b2),ylim=c(0,70),xlim=c(0,0.4),type="l")
# 

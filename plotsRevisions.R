# Graphs copula

library(copula)
library(HAC)

source("/media/baeyc/Seagate Backup Plus Drive/CEC/Code/es2cropBB/pollination/calibration/floralCoverage/helpers.R")

mu1<-0.1
mu2<-0.05
s1 <- 50
s2 <- 30

tauB<-0.2
tauW<-0.6

l1 <- paramBetaFromModeSampSize(mu1,s1)
l2 <- paramBetaFromModeSampSize(mu2,s2)

thetaW<-copJoe@iTau(tauW)
thetaB<-copJoe@iTau(tauB)

# we define the nested copula, type=1 -> Gumbel, 3 -> Clayton, 5 -> Frank and 7 -> Joe
cop<-hac(type=7,tree=list(list("y1","y2",thetaW),list("y3","y4",thetaW),thetaW))

vecSimu <- rHAC(5000,cop)
vecBeta <- 0*vecSimu
vecBeta[,1] <- qbeta(vecSimu[,1],l1[[1]],l1[[2]])
vecBeta[,2] <- qbeta(vecSimu[,2],l1[[1]],l1[[2]])
vecBeta[,3] <- qbeta(vecSimu[,3],l2[[1]],l2[[2]])
vecBeta[,4] <- qbeta(vecSimu[,4],l2[[1]],l2[[2]])

require(GGally)
pairs(vecBeta,pch=20)

d <- as.data.frame(vecBeta)
p <- ggpairs(d, lower = list(continuous = wrap("points", alpha = 0.3)), upper = list(continuous = 'blank'), diag = list(continuous = 'blank'))

ggsave("~/Code/es2crop/pollination/calibration/floralCoverage/pairsJoe.pdf", p, height = 5)


# Kendall's tau real data
source("/media/baeyc/Seagate Backup Plus Drive/CEC/Code/es2cropBB/pollination/calibration/floralCoverage/loadObsData.R")
require(VGAM)

kendall.tau(osr[,1],osr[,2])

# field edges
# within period
kendall.tau(fieldEdge[,1],fieldEdge[,2])
kendall.tau(fieldEdge[,3],fieldEdge[,4])
# between period
kendall.tau(fieldEdge[,1],fieldEdge[,3])
kendall.tau(fieldEdge[,1],fieldEdge[,4])
kendall.tau(fieldEdge[,2],fieldEdge[,3])
kendall.tau(fieldEdge[,2],fieldEdge[,4])

# field edges with flowers
# within period
kendall.tau(fieldEdgeFl[,1],fieldEdgeFl[,2])
kendall.tau(fieldEdgeFl[,3],fieldEdgeFl[,4])
# between period
kendall.tau(fieldEdgeFl[,1],fieldEdgeFl[,3])
kendall.tau(fieldEdgeFl[,1],fieldEdgeFl[,4])
kendall.tau(fieldEdgeFl[,2],fieldEdgeFl[,3])
kendall.tau(fieldEdgeFl[,2],fieldEdgeFl[,4])

# pasture
# within period
kendall.tau(pasture[,1],pasture[,2])
kendall.tau(pasture[,3],pasture[,4])
# between period
kendall.tau(pasture[,1],pasture[,3])
kendall.tau(pasture[,1],pasture[,4])
kendall.tau(pasture[,2],pasture[,3])
kendall.tau(pasture[,2],pasture[,4])



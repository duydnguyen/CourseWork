rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()
library(ape)
#### This is STAT877 HW_02 ####

#### Initialize Header File ####
FilePath <- '~/github/STAT877_Spring2014/HW_02'
################################

################ Input data #################
setwd(FilePath)
source("phylogeny.r")
Filename <- "HIVSIVgag.nex"
xlist = read.nexus.data(Filename)
xbin = as.DNAbin(xlist)


############### Question 1: Maximum Likelihood ############
# ordered: A, C, G, T
p = c(0.2,0.3,0.4,0.1) 
kappaR=3 # r_AG
kappaY=4 # r_CT
Q <- TN93(p=p,kappaR=kappaR,kappaY=kappaY)
t <- 0.03
n12 <- table(xlist[[1]],xlist[[2]])
P.t <- matrixExp(Q,t)
log.lh = 0
for (i in 1:4) {
    for (j in 1:4) {
        log.lh = log.lh + n12[i,j]*( log(p[i]) + log(P.t[i,j]))          
    }
}

Eval.Loglh <- function(p, kappaR, kappaY, t, seq1, seq2){
  Q <- TN93(p=p,kappaR=kappaR,kappaY=kappaY)  
  n12 <- table(seq1,seq2)
  P.t <- matrixExp(Q,t)
  log.lh = 0
  for (i in 1:4) {
    for (j in 1:4) {
      log.lh = log.lh + n12[i,j]*( log(p[i]) + log(P.t[i,j]))          
    }
  }
  return(log.lh)
}

try = Eval.Loglh(c(0.2,0.3,0.4,0.1), 3, 4, 0.03, xlist[[1]], xlist[[2]])








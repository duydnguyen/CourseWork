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
####### (a) ######
# ordered: A, C, G, T
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

#try = Eval.Loglh(c(0.2,0.3,0.4,0.1), 3, 4, 0.03, xlist[[1]], xlist[[2]])
## Compute p given two sequences. The base.freq() does not work since it varies every time I run
Eval.pi <- function(seq1, seq2){
  p=c(rep(0,4))
  total = length(which(seq1 !='-')) + length(which(seq2 !='-'))
  p[1] <- ( length(which(seq1 =='a')) + length(which(seq2 =='a')) )/ total
  p[2] <- ( length(which(seq1 =='c')) + length(which(seq2 =='c')) )/ total
  p[3] <- ( length(which(seq1 =='g')) + length(which(seq2 =='g')) )/ total
  p[4] <- ( length(which(seq1 =='t')) + length(which(seq2 =='t')) )/ total
  return(p)  
}
## BALT86 (HIVstrain) and sm83 (SIV strain) with kappaR=3, kappaY=4
p = Eval.pi(xlist[[3]], xlist[[1]])
# time starts backward from 1983?
time = seq(0.5,50,0.5)
log.lh = c(rep(0, length(time)))
for (i in 1:length(time)){
    log.lh[i] = Eval.Loglh(p, 3, 4, time[i], xlist[[3]], xlist[[1]])
}
plot(log.lh ~ time, xlab="function of distance", ylab="log-likelihood", type="l")
####### (b) ######
## Compute pi for entire data
seq1 = xlist[[1]]
seq2 = xlist[[2]]
seq3 = xlist[[3]]
seq4 = xlist[[4]]
p=c(rep(0,4))
total = length(which(seq1 !='-')) + length(which(seq2 !='-')) +
          length(which(seq3 !='-')) + length(which(seq4 !='-'))
p[1] <- ( length(which(seq1 =='a')) + length(which(seq2 =='a')) 
          + length(which(seq3 =='a')) + length(which(seq4 =='a')) ) / total
p[2] <- ( length(which(seq1 =='c')) + length(which(seq2 =='c')) 
          + length(which(seq3 =='c')) + length(which(seq4 =='c')) ) / total
p[3] <- ( length(which(seq1 =='g')) + length(which(seq2 =='g')) 
          + length(which(seq3 =='g')) + length(which(seq4 =='g')) ) / total
p[4] <- ( length(which(seq1 =='t')) + length(which(seq2 =='t')) 
          + length(which(seq3 =='t')) + length(which(seq4 =='t')) ) / total


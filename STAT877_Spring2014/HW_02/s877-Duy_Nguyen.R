rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()
library(ape)
#### STAT877 HW_02 ####

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
  d = dim(n12)
  log.lh = 0
# No missing values  
  if ((d[1]==4)&(d[2]==4))
    for (i in 1:4) {
      for (j in 1:4) {
        log.lh = log.lh + n12[i,j]*( log(p[i]) + log(P.t[i,j]))          
      }
    }
# Missing values
if ((d[1]==5)&(d[2]==4))
{
  for (j in 1:4)
  {
    log.lh = log.lh + n12[1,j]*(log(p[j]))
  }
  for (i in 2:5)
    for (j in 1:4)
    {
      log.lh = log.lh + n12[i,j]*(log(p[i-1])+log(P.t[i-1,j]))
    }
}

if ((d[1]==4)&(d[2]==5))
{
  for (i in 1:4)
  {
    log.lh = log.lh + n12[i,1]*(log(p[i]))
  }
  for (i in 1:4)
    for (j in 2:5)
    {
      log.lh = log.lh + n12[i,j]*(log(p[i])+log(P.t[i,j-1]))
    }
}
return(log.lh)
}

## Compute p given the entire data. 
## The base.freq() does not work since it varies every time I run
Eval.pi <- function(){
    # Declare my vector pi p()
    p=c(rep(0,4))
    df <- c(xlist[[1]],xlist[[2]],xlist[[3]],xlist[[4]])
    index = which(df=="-")
    if (length(index)>0) df = df[-index]
    p[1]=sum(df=="a")/length(df)
    p[2]=sum(df=="c")/length(df)
    p[3]=sum(df=="g")/length(df)
    p[4]=sum(df=="t")/length(df)
    return(p)  
}
# Find p 
p = Eval.pi()
## BALT86 (HIVstrain) and sm83 (SIV strain) with kappaR=3, kappaY=4

time = seq(0.01,1,0.001)
log.lh = c(rep(0, length(time)))
for (i in 1:length(time)){
    log.lh[i] = Eval.Loglh(p, 3, 4, time[i], xlist[[3]], xlist[[1]])
}
#Plot likelihood function for distance d(sm3, BALT86)
plot(log.lh ~ time,xlab="t",ylab="Likelihood",main="Likelihood function for d(sm83,BALT86)",
     type='l',cex.lab=1.5, cex.axis=1.5,cex.main=1.5)

####### (b) ######
## function to find maximum likelihood distance given each pair of sequences
Eval.max <- function(seq1, seq2){
  fn.b = function(t){return(Eval.Loglh(p,2,3,t, seq1, seq2))}
  MLE = optim(0.1,fn.b,control=list(fnscale=-1),method="Brent",lower=0.001,upper=1)
  return(MLE)
} 

d12 <- Eval.max(xlist[[1]],xlist[[2]])
d13 <- Eval.max(xlist[[1]],xlist[[3]])
d14 <- Eval.max(xlist[[1]],xlist[[4]])
d23 <- Eval.max(xlist[[2]],xlist[[3]])
d24 <- Eval.max(xlist[[2]],xlist[[4]])
d34 <- Eval.max(xlist[[3]],xlist[[4]])

# To get the maximum likelihood distance, use d12$par

####### (c) ######
## function in R that will optimize the likelihood: parameters = distance, kappaR, kappaY
## x[1]=kappaR, x[2]=kappaY, x[3]= distance
Eval.max.global <- function(seq1,seq2){
  fn.c <- function(x){return(Eval.Loglh(p, x[1], x[2], x[3], seq1, seq2 ))}
  MLE.c <- optim(c(2,2.8,0.1),fn.c,control=list(fnscale=-1))
  return(MLE.c)    
} 

seq1 <- xlist[[1]]
seq2 <- xlist[[3]]

d12.global <- Eval.max.global(xlist[[1]],xlist[[2]])
# MLEs: kappaR, kappaY, distance(time)
d12.global$par

d13.global <- Eval.max.global(xlist[[1]],xlist[[3]])
d14.global <- Eval.max.global(xlist[[1]],xlist[[4]])
d23.global <- Eval.max.global(xlist[[2]],xlist[[3]])
d24.global <- Eval.max.global(xlist[[2]],xlist[[4]])
d34.global <- Eval.max.global(xlist[[3]],xlist[[4]])

## Plots of likelihood: I plot d12. Other cases are similar
MLE <- d12.global$par
y = c()
x = seq(0.01,1,0.001)
for (i in 1:length(x))
{
  y[i] = Eval.Loglh(p,MLE[1],MLE[2],x[i],xlist[[1]],xlist[[2]])
}
plot(y~x, type='l')
############### Question 2: MCMC ############
## Proposal distribution q =Exp(t)
kappaR = 2
kappaY = 3
seq1 <- xlist[[1]]
seq2 <- xlist[[3]]
p = Eval.pi(seq1, seq2)
# target distribution:
f.pr <- function(t)
{
  return(10*exp(-10*t + Eval.Loglh(p, kappaR,kappaY,t,seq1,seq2)))
} 
## Store time distance into vector t 
# N = number of samples
Sim <-function(N){
  t= c()
  t[1]=0.2
  for(i in 2:N)
  {
    tnew =rexp(1,rate = 4)  ## can be fixed
  
   b = f.pr(tnew)*dexp(t[i-1],rate = tnew)/(f.pr(t[(i-1)])*dexp(tnew,rate = t[i-1]))
   p.accept = min(1,b)
   index = rbinom(1, size = 1, prob = p.accept)
   if (index==1) t[i]=tnew
   else t[i]=t[i-1]
  }
  return(t)
}

t <- Sim(100000)
t
## Plots

plot(t[99001:10000],type='l',main=' Trace plot for MCMC method',ylab="t",cex.lab=1.5,
     cex.axis=1.5,cex.main=1.5)
hist(t[99001:10000],main='Posterior distribution',freq=FALSE, xlab = "t",cex.lab=1.5,
     cex.axis=1.5,cex.main=1.5)


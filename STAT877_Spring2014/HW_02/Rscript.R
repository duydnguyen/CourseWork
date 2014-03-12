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
n12 = table(xlist[[1]],xlist[[2]])

############### Question 1: Maximum Likelihood ############












rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

#### Initialize Header File ####
FilePath <- '~/github/STAT877_Spring2014/qtlproject'
Filename.Header <- paste('~/RScripts/HeaderFile_HW.R', sep='')
source(Filename.Header)
source(paste(FilePath, 'fn_Library.R', sep=''))
################################

################ Input data #################
setwd(FilePath)
Filename <- 'hw.dat'
Data <- read.table(file=Filename, header=FALSE)

################# Marker Regression #########



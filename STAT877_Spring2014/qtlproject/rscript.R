rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

#### This is HW? ####

#### Initialize Header File ####
FilePath <- '~/github/STAT850_Spring2014/'
Filename.Header <- paste('~/RScripts/HeaderFile_HW.R', sep='')
source(Filename.Header)
source(paste(FilePath, 'fn_Library.R', sep=''))
################################

################ Input data #################
setwd(FilePath)
Filename <- 'hw.dat'
Data <- read.table(file=Filename, header=FALSE)


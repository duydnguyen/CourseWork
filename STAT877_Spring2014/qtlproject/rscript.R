rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

#### Initialize Header File ####
FilePath <- '~/github/STAT877_Spring2014/qtlproject'
#Filename.Header <- paste('~/RScripts/HeaderFile_HW.R', sep='')
#source(Filename.Header)
#source(paste(FilePath, 'fn_Library.R', sep=''))
################################

################ Input data #################
setwd(FilePath)
#Filename <- 'hw.dat'
#Data <- read.table(file=Filename, header=FALSE)
library(qtl)
library(flexmix)
data(hyper)
str(hyper)
summary(hyper)
################# Marker Regression #########
out.mr <- scanone(hyper, method="mr")
out.mr[out.mr$chr == 4,]
# focus on chr 4, marker D4Mit214. Chr 4 has 20 markers. Codes: 1=homozygotes, 2=heterozygotes
hyper$geno[[4]]$map[1]
hyper$geno[[4]]
str(hyper$geno[[4]])
hyper$pheno[1:5,]
hyper$geno[[4]]$data[1:3,]
str(hyper$geno[[4]]$data)
### Extract yn=phenotype, x=genotypes on chr 4 at marker D4Mit214 (=col 6 in data). No missing data
data.chr4 <-hyper$geno[[4]]$data
yn <- hyper$pheno$bp
x <- as.factor(data.chr4[,6])
boxplot(yn ~ x)
## t-test
#anova( lm(yn ~ x) )
#t.test(yn ~ data.chr4[,6], alternative=c("two.sided"), var.equal=TRUE)
model.chr4 <- aov(yn ~ x)
summary(model.chr4) #F=33.43 -> LOD = 6.8648387 same as in out.mr
## Perform LCR: yn ~ x. Remark: coefficients in each class vary every time run E-M ???
model.LCR <- flexmix(yn ~ x, k =2 )
model.LCR
parameters(model.LCR, component = 1)
parameters(model.LCR, component = 2)
# labels for yn
class1 <- which(model.LCR@cluster == 1)
class2 <- which(model.LCR@cluster == 2)
# Test for significance
rmodel.LCR <- refit(model.LCR)
summary(rmodel.LCR)



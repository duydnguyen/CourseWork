T_real_mdd <- read.table("test_scores_real_mdd", quote="\"")
T_false_mdd <- read.table("test_scores_false_mdd", quote="\"")
T_real <- read.table("test_scores_real", quote="\"")
T_false <- read.table("test_scores_false", quote="\"")

df <- data.frame(T_real, T_false, T_real_mdd, T_false_mdd)
colnames(df) <- c('realPWM', 'falsePWM', 'realMDD', 'falseMDD')
View(df)

#hist(rbind(df$realPWM, df$falsePWM))
#hist(rbind(df$realMDD, df$falseMDD))
library(ROCR)
### Plot PC curve
real = matrix(df$realMDD, nrow = 1000, ncol =1)
false = matrix(df$falseMDD, nrow = 1000, ncol =1)
scores = rbind(real, false)
min = min(scores)
max = max(scores)
T_seq = seq(min, max, length.out = 100)

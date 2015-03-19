T_real_mdd <- read.table("test_scores_real_mdd", quote="\"")
T_false_mdd <- read.table("test_scores_false_mdd", quote="\"")
T_real <- read.table("test_scores_real", quote="\"")
T_false <- read.table("test_scores_false", quote="\"")

df <- data.frame(T_real, T_false, T_real_mdd, T_false_mdd)
colnames(df) <- c('realPWM', 'falsePWM', 'realMDD', 'falseMDD')


getLabel <- function(scores, cutoff) {
  length = dim(scores)[1]
  lab_pred <- rep(NA, length )
  for (i in 1:length) {
    if (scores[i] > cutoff) lab_pred[i] <- 1
    else lab_pred[i] <- 0
   }  
  return(lab_pred)
}


### Plot PC curve: labels (real = 1, false = 0 )
evalPC <- function(real, false) {
  real = matrix(real, nrow = 1000, ncol =1)
  false = matrix(false, nrow = 1000, ncol =1)
  scores = rbind(real, false)
  min = min(scores)
  max = max(scores)
  lenT = 100
  T_seq = seq(min, max, length.out = lenT)
  Precision = rep(0, lenT)
  Recall = rep(0, lenT)
  for (i in 1:lenT) {
    cutoff = T_seq[i]
    lab_pred <- getLabel(scores, cutoff)
    ## Compute TP
    TP = length(which(lab_pred[1:1000] == 1))
    FN = length(which(lab_pred[1:1000] == 0))
    FP = length(which(lab_pred[1001:2000] == 1))
    Precision[i] = TP / (TP + FP)
    Recall[i] = TP / (TP + FN)  
  }
  PC <- data.frame(Precision, Recall)
  colnames(PC) <- c('Precision', 'Recall')
  return(PC)
}

PCmdd <- evalPC(df$realMDD, df$falseMDD)
PCpwd <- evalPC(df$realPWM, df$falsePWM)
# Plots
plot(x = PCpwd$Recall, y = PCpwd$Precision, type = 'l', col = 1,
     xlab = 'Recall', ylab = 'Precision', main = 'Precision-Recall Curve')
points(x = PCmdd$Recall, y = PCmdd$Precision, type = 'l', col = 2 )
legend(0.2,0.6, c('PWD', 'MDD'), lty = c(1,1), lwd = c(2.5,2.5), col = c(1,2))



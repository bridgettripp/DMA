## This code identifies the ideal strength cutoff for learned edges
## Adapted from Scutari M and Nagarajan R (2013) 'Identifying significant
   ## edges in graphical models of molecular networks'
## input: column vector of strength values
## output: strength cutoff for identifying true edges

library(plyr)

confidenceValue <- function(x){
  str.sort <- as.data.frame(sort(round(x, 3)))
  str.sort <- as.data.frame(as.data.frame(split(str.sort,1:2))[,1])
  colnames(str.sort) <- 'strength'
  str.count <- count(str.sort, 'strength')
  int.total = dim(str.sort)[1]
  L <- data.frame()
  for (e in seq(dim(str.count)[1])){
    str.count[e,3] <- str.count[e+1,1]
    str.count[e,4] <- str.count[e,3] - str.count[e,1]
    str.count[e,5] <- (sum(str.count[1:e,2])/int.total)
  }
  str.count[dim(str.count)[1], 3:4] <- c(1, 0)
  for (t in seq(0,1,by = (1/(int.total*4)))){ ###10000
    L <- rbind.data.frame(L,t)} 
  v <- vector()
  t.hat <- vector()
  for (t in seq(dim(L)[1])){
    v <- apply(str.count,1,function(y) (abs(y[5]-L[t,]) * y[4]))
    v <- as.data.frame(v)
    t.h <- sum(v[1:dim(v)[1],])
    t.hat <- rbind.data.frame(t.hat, t.h) 
  }
  zebra <- cbind(L, t.hat)
  colnames(zebra) <- c("T","L")
  colnames(str.count) <- c("value1", "freq", "value2", "diff(v2-v1)","PDF")
  min.t.hat <- zebra[which.min(zebra[,2]),1]
  min.l <- (zebra[which.min(zebra[,2]),2]*10)
  str.count[which(min.t.hat < str.count[,5]),1][1]
}
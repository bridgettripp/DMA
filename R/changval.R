
changval <- function(ct, cs, p){
  ct <- rbind.data.frame(cs[i,], ct); ct[1,p] <- NA
  ct2 <- impute::impute.knn(as.matrix(ct), k = round(sqrt(dim(ct)[1])))[["data"]][1,] # try larger k
  #ct <- ct[-1]
  return(ct2)
}


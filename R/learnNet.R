learnStructure <- function(input_data){
  set.seed(3080)
  source("../../network-optimization-tools/confidenceValue.R")
  cl = parallel::makeCluster(4)
  z <- round(.85*(dim(input_data[1])))[1]
  strength <- bnlearn::boot.strength(input_data, R = 1000, 
                                   m = z, algorithm = 'mmhc',
                                   cpdag = TRUE, cluster = cl) 

  parallel::stopCluster(cl)
  #strg <- (confidenceValue(strength$strength)*.75)
  #avg.boot <- bnlearn::averaged.network(strength, threshold = strg)
  return(strength)
}


setwd("~/Dropbox/dissertation/Aim2/network_alteration")

rm(list=ls())

source("calcSHD.R")
source("changval.R")
source("learnNet.R")
source("http://www.sthda.com/upload/rquery_cormat.r")
source("~/Dropbox/Codes/OBaNK/confidenceValue.R")

library(bnlearn)
library(tidyr)
library(plyr)

#' load files

data0 <- read.csv("~/Dropbox/dissertation/Aim2/datasets/pregnancy_cnt.csv", check.names = FALSE)
data0 <- data0 %>%
  dplyr::select(-Symbol)
case0 <- read.csv("~/Dropbox/dissertation/Aim2/datasets/pregnancy_tri3.csv", check.names = FALSE)
case0 <- case0 %>%
  dplyr::select(-Symbol)
prefx <- "sandbox"
mol_type <- read.csv("~/Dropbox/dissertation/Aim2/datasets/pregnancy_types.csv", check.names = FALSE)

#' jitter
set.seed(0328)
library(FNN)
## control data
data1 <- data0 
#x = length(data0) #118
x = 118
n <- dim(data0)[1]
for (i in 1:dim(data0)[2]){
  m1 <- median(unlist(data0[,i])) 
  sd1 <- diff(range(data0[,i]))/4
  ri <- range(data0[,i])
  for (j in 1:x) {
    data1[n+j, i] <- rnorm(1, m1, sd1)
  }}
gh <- KL.divergence(as.numeric(unlist(data0)), as.numeric(unlist(data1)), k = 50, algorithm="kd_tree")
kld <- round(mean(gh[5:50]), 3)
data0 <- data1

## case data
set.seed(0328)
data1 <- case0 
#x = length(data0)
x = 118
n <- dim(case0)[1]
for (i in 1:dim(case0)[2]){
  m1 <- median(unlist(case0[,i])) 
  sd1 <- diff(range(case0[,i]))/4
  ri <- range(case0[,i])
  for (j in 1:x) {
    data1[n+j, i] <- rnorm(1, m1, sd1)
  }}
case0 <- data1


######################
#' learn structure 
#' base/control network
strength <- learnStructure(data0)
strg <- (confidenceValue(strength$strength)*.75) 
avg.boot <- bnlearn::averaged.network(strength, threshold = strg)


#' active/case network
#' Network that is actively altered 
strength2 <- learnStructure(case0)
strg2 <- (confidenceValue(strength2$strength)*.75)
avg.boot2 <- bnlearn::averaged.network(strength2, threshold = strg2)


#' plot original base and active networks
pdf(paste0(prefx, "_intial.pdf"))
par(mfrow = c(1,1), font = 2, pty="m")
bnlearn::graphviz.compare(bnlearn::skeleton(avg.boot), 
                          bnlearn::skeleton(avg.boot2), diff = "none", 
                          main = c("Control Network", "Case Network"), layout = "fdp", shape = "rectangle")
dev.off()

par(mfrow=c(1,2), font = 2, pty="m")
bnlearn::graphviz.compare(bnlearn::skeleton(avg.boot), 
                          bnlearn::skeleton(avg.boot2), diff = "none", main = c("Control Network", "Case Network"), layout = "fdp", shape = "rectangle")

#####################################

#' Isolate nodes with 1 or greater links.
#' Removing isolated nodes for centrality calculation
nodes <- unique(unlist(arcs(avg.boot)))
#averaged <- avg.boot

relnodes = nodes(avg.boot)[sapply(nodes, degree, object = avg.boot) > 0]
relnodes = relnodes[!is.na(relnodes)]

str2 = strength[(strength$from %in% relnodes) & (strength$to %in% relnodes) & 
                  (strength$strength > strg), ]

avg.boot0 <- bnlearn::averaged.network(str2, threshold = strg)

source("testCentrality.R")
#centrality <- calc_cent_weighted(avg.boot, mol_type)
centrality <- calc_cent(avg.boot)

sink(paste0(prefx,"-centrality.txt"), append = TRUE)
print(centrality)
sink()


org_auc <- calc_auc(avg.boot, strength2)
org_frobnorm <- frob_norm(strength, strength2)

frobnorm <- org_frobnorm

plot_fbn(avg.boot, avg.boot2, frobnorm, "Original Case")

current_frobnorm <- org_frobnorm
current_auc <- org_auc

original_case0 <- case0
current_case0 <- case0

high_auc <- org_auc
low.fbn <- org_frobnorm

nodeorder <- vector()

fbn <- tibble()
fbn <- rbind(centrality[1], fbn); colnames(fbn) <- "node"

cent_list <- as.data.frame(centrality[1]); colnames(cent_list) <- "node"

altered_node <- "NA"

#' Account for networks with fewer nodes than 15
if(dim(centrality)[1] < 15){
  end_num <- dim(centrality)[1]
} else {
  end_num <- 15
}

set.seed(03)
for (h in 1:end_num){
  #' remove temp files to free up previously used files
  tmp_dir <- tempdir()
  file.remove(list.files(tmp_dir, full.names = T,
                         pattern = "^file"))
  print(c("start:", h))
  for (j in 1:end_num){
    print(j)
    modified_case_net <- data.frame(matrix(ncol = length(current_case0), nrow = 0))
    colnames(modified_case_net) <- colnames(current_case0)
    for (i in 1:dim(current_case0)[1]){
      #' Find column name that matches next centrality node.
      #' Alter signal values
      modified_case_net <- rbind.data.frame(modified_case_net, changval(data0, current_case0,
                                              which(cent_list[j,1]==colnames(current_case0))))
      colnames(modified_case_net) <- colnames(current_case0)
    }
    
    #' Learn structure with new imputed values
      print("learn structure with new values")
      print(j)
      strength2.5 <- learnStructure(modified_case_net)
      strg <- (confidenceValue(strength2.5$strength)*.75)
      avg.boot2.5 <- bnlearn::averaged.network(strength2.5, threshold = strg)
      
    #' calculate frobenius norm of the difference between strengths
      print("calc frobenius")
      frobnorm <- frob_norm(strength, strength2.5)
 
    #' calculate the auc using original and new network
      print("calc auc")
      ca1 <- calc_auc(avg.boot, strength2.5)

      if(frobnorm <= (low.fbn*1.03)){ #Slight increase is to limit ending prematurely in local minimum
    #' If new fbn is lower than previous, new network replaces previous
        print(cent_list[j,1])
        low.fbn <- frobnorm
        print(low.fbn)
        best_case <- modified_case_net
        altered_node <- cent_list[j,1]
        print(altered_node)
        avg.bootBest <- avg.boot2.5
    }
  }
  
  fbn[which(fbn$node==cent_list[j,1]), 2:3] <- c(frobnorm, ca1)
  colnames(fbn) <- c("node", "fbn", "auc")
  
  pdf(paste0(prefx,"-","node", altered_node, ".pdf"))
  par(mfrow = c(1,1), font = 2, pty="m")
  plot_fbn3(avg.boot, avg.boot2,avg.bootBest, low.fbn, altered_node)
  par(mfrow = c(1,3), font = 2, pty="m")
  plot_fbn3(avg.boot, avg.boot2,avg.bootBest, low.fbn, altered_node)

  dev.off()

  cent_list <- as.data.frame(cent_list[-which(cent_list==altered_node),1])
  

  pdf(paste0(prefx,"-final.pdf"), family = "Times")
  par(mfrow = c(1,3), font = 2, pty="m")
  plot_fbn(avg.boot, avg.bootBest, low.fbn, altered_node)
  plot_fbn3(avg.boot, avg.boot2, avg.bootBest, low.fbn, altered_node)
  dev.off()

  current_case0 <- best_case
  print(paste("end:", h))
  
  print(paste("added", cent_list[j,1]))

  if(altered_node %in% nodeorder == FALSE){
    end_num <- (end_num - 1)
    nodeorder <- cbind(nodeorder, altered_node)

  } else
    {
    print("done")
    break }
}

#' remove temp files to free up previously used files
tmp_dir <- tempdir()
file.remove(list.files(tmp_dir, full.names = T, 
                       pattern = "^file"))

write.csv(nodeorder, "node_order.csv")
save.image(file = paste0(prefx, "-DMA_final.rda"))


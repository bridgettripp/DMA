
calc_shd <- function(true_net, learned_net){
  shd_score <- bnlearn::shd(bnlearn::skeleton(true_net), bnlearn::skeleton(learned_net), debug = TRUE)
  return(shd_score)
}

plot_shd <- function(true_net, learned_net, cshd){
  par(mfrow = c(1,2), font = 2, pty="m")
  bnlearn::graphviz.compare(bnlearn::skeleton(true_net),
                          bnlearn::skeleton(learned_net), diff = "none",
                          main = c("Control Network", paste("Case Network", "SHD:", cshd)), layout = "fdp", shape = "rectangle")}
                     

frob_norm <- function(cnt, case){
  cnt$combo <- paste(cnt$from, cnt$to, sep = "_")
  case$combo <- paste(case$from, case$to, sep = "_")
  c2 <- dplyr::inner_join(cnt, case, by = "combo"); n <- dim(c2)[1]
  #fbm <- round(sqrt(sum(as.matrix((c2$strength.x-c2$strength.y)^2))), digits = 4)
  #             #/n, digits = 6)
  fbm <- norm(as.matrix(c2$strength.x-c2$strength.y), "F")
  #fbm <- abs(norm(as.matrix(c2$strength.x), "F")-norm(as.matrix(c2$strength.y), "F"))
 
}


plot_fbn <- function(true_net, learned_net, fbn, altered_node){
  par(mfrow = c(1,2), font = 2, pty="m")
  bnlearn::graphviz.compare(bnlearn::skeleton(true_net),
                            bnlearn::skeleton(learned_net), diff = "none",
                            main = c("Control Network", paste(altered_node, round(fbn, 3))), layout = "fdp", shape = "rectangle")}

plot_fbn3 <- function(true_net, case_net, learned_net, fbn, altered_node){
  par(mfrow = c(1,3), font = 2, pty="m")
  bnlearn::graphviz.compare(bnlearn::skeleton(true_net), bnlearn::skeleton(case_net),
                            bnlearn::skeleton(learned_net), diff = "none",
                            main = c("Control Network", "Original Case", paste(altered_node, round(fbn, 3))), layout = "fdp", shape = "rectangle")}


calc_auc <- function(true.net, learned.strength){
  true.arcs <- true.net$arcs #true.net is a bnlearn object
  learned.arcs <- as.data.frame(cbind(learned.strength$from, learned.strength$to, 
                                      learned.strength$strength)) 
  colnames(learned.arcs) <- c("from","to", "streng")
  strength1 <- learned.strength
  
  for (i in 1:dim(learned.strength)[1]){
    if (length(which(learned.strength$from[i]==true.arcs[,1] & learned.strength$to[i]
                     ==true.arcs[,2])) == 1) {
      learned.strength[i,5] <- 1}
    else {
      learned.strength[i,5] <- 0
    }}
  
  predictions <- as.data.frame(learned.strength$strength)
  labels <- as.matrix(learned.strength[,5])
  pred <- ROCR::prediction(predictions, labels)
  perf <- ROCR::performance(pred,'tpr','fpr')
  auc.perf <- ROCR::performance(pred, "auc")
  auc1 <- as.numeric(auc.perf@y.values)
  return(auc1)
}
calc_cent <- function(netw){
  b1 <- bnlearn::as.igraph(bnlearn::skeleton(netw))
  rmedges <- names(which(igraph::degree(b1)==0))
  b1 <- igraph::delete.vertices(b1,c(rmedges))
  components <- igraph::clusters(b1, mode="weak")
  biggest_cluster_id <- which.max(components$csize)
  node_ids <- igraph::V(b1)[components$membership == biggest_cluster_id]
  b1 <- igraph::induced_subgraph(b1, node_ids) #induced subgraph
  #
  p <- CINNA::pca_centralities( CINNA::calculate_centralities(b1, include = c(
                                                                              "Closeness Centrality (Freeman)",
                                                                              "Degree Centrality",
                                                                              "eigenvector centralities",
                                                                              "Markov Centrality",
                                                                              "Laplacian Centrality"
                                                                              )), scale = FALSE )[["data"]]
  bcent <- gsub(".", " ", as.matrix(p[which(p$contrib == max(p$contrib)),1]), fixed = TRUE)
  cent_table <- as.data.frame(CINNA::calculate_centralities(b1, include = bcent))
  cent_table <- cbind.data.frame(bnlearn::nodes(as.bn(b1)), cent_table)
  cent_table <- cent_table[order(-cent_table[,2], decreasing = FALSE), ]
  return(cent_table)
}


calc_cent_weighted <- function(netw, mol_type){
  library(data.table)
  moltype <- mol_type; moltype <- as.data.frame(mol_type) 
  b1 <- bnlearn::as.igraph(bnlearn::skeleton(netw))
  rmedges <- names(which(igraph::degree(b1)==0))
  b1 <- igraph::delete.vertices(b1,c(rmedges))
  components <- igraph::clusters(b1, mode="weak")
  biggest_cluster_id <- which.max(components$csize)
  node_ids <- igraph::V(b1)[components$membership == biggest_cluster_id]
  b1 <- igraph::induced_subgraph(b1, node_ids) #induced subgraph
  b2 <- igraph::as_data_frame(b1) 
  for (m in 1:nrow(b2)){
    b2$to[m] <- moltype[moltype==b2$to[m],2]
    b2$main[m] <- moltype[moltype==b2$from[m],2]
    b2$match[m] <-ifelse(b2$main[m]==b2$to[m], 0, 1)
    }
  setDT(b2)
  b2 <- b2[ ,list(sum=sum(match)), by=from]; colnames(b2) <- c("node", "ttl") 
  
  p <- CINNA::pca_centralities( CINNA::calculate_centralities(b1, include = c(
    "Closeness Centrality (Freeman)",
    "Degree Centrality",
    "eigenvector centralities",
    "Markov Centrality",
    "Laplacian Centrality"
  )), scale = FALSE )[["data"]]
  
  bcent <- gsub(".", " ", as.matrix(p[which(p$contrib == max(p$contrib)),1]), fixed = TRUE)
  cent_table <- as.data.frame(CINNA::calculate_centralities(b1, include = bcent))
  cent_table <- cbind.data.frame(bnlearn::nodes(as.bn(b1)), cent_table)
  cent_table <- cent_table[order(-cent_table[,2], decreasing = FALSE), ]
  colnames(cent_table) <- c("node", "centrality")
  sd1 <- as.numeric(sd(cent_table$centrality)/2)
  cent_table <- cent_table %>%
    dplyr::inner_join(b2) %>%
    mutate(factor = ((ttl * sd1) + centrality)) %>%
    mutate(centrality = dplyr::if_else(factor >= centrality, factor, centrality))
  cent_table <- cent_table[,-(3:4)]
  cent_table <- cent_table[order(cent_table$centrality, decreasing = TRUE), ]
  return(cent_table)
}



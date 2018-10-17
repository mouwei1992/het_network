pairwise_dist <- function(positions){
  n_nodes = ncol(positions)
  dist_matrix = matrix(0, ncol = n_nodes, nrow = n_nodes)
  
  for(node_1 in 1:n_nodes){
    for(node_2 in 1:n_nodes){
      if(node_1 >= node_2)  dist_matrix[node_1, node_2] = 0
      else  dist_matrix[node_1, node_2] = norm(as.matrix(positions[, node_1] - positions[, node_2]), type = "F")
    }
  }
  
  return(dist_matrix)
}

dist_ratio <- function(pos_1, pos_2){
  pairwise_dist_1 = pairwise_dist(pos_1)
  pairwise_dist_2 = pairwise_dist(pos_2)
  
  dist_ratio = pairwise_dist_1 / pairwise_dist_2
  
  return(dist_ratio[pairwise_dist_2 != 0])
}

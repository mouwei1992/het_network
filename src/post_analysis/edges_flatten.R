edges_flatten <- function(adj_mat, latent_pos.x, latent_pos.y){
  n_nodes <- ncol(adj_mat)
  flattened <- matrix(ncol = 4, nrow = 0)
  
  row_count = 1
  for(i in 1:n_nodes){
    for(j in 1:n_nodes){
      if(i==j) next
      if(adj_mat[i,j] == 1){
        flattened = rbind(flattened, 
                          c(latent_pos.x[i],latent_pos.y[i],
                            latent_pos.x[j],latent_pos.y[j]))
      }
    }
  }
  
  flattened <- data.frame(flattened)
  names(flattened) <- c("x","y","xend","yend")
  return(flattened)
}

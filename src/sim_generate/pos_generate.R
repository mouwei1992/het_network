setwd("~/Dropbox/het_network_project/version6/src/sim_generate")

# generate the true position
position_gen <- function(n_clusters, n_nodes, space_dim){
  true_pos = 
    matrix(rnorm(space_dim * n_nodes), nrow = space_dim, ncol = n_nodes)
  clusters_centers = matrix(0, nrow = space_dim, ncol = n_cluseters)
  
  
  return true_pos;
}

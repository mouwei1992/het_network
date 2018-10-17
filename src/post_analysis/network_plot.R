library(ggplot2)

# initialization
positions2.data.frame <- 
  function(latent_positions, nodes.type = NULL){
  return(data.frame(
    x = latent_positions[1,],
    y = latent_positions[2,],
    nodes.type = as.factor(nodes.type)
    ))
}

init_canvas <- function(latent_positions, nodes.type){
  ggplot(positions2.data.frame(latent_positions, nodes.type)) + 
    theme_minimal()
}

# plotting nodes..
plot_points <- function(ggplot.canvas){
  (ggplot.canvas + 
    geom_point(aes(
      x = x,
      y = y,
      color = nodes.type
      ),
      size = 1
      ))
}

# 


# for undirected graph
# this function cannot plot too many edges at one time
# in that case, run it multiple times or use loops

plot_edges <- function(ggplot.canvas, adj_mat){
  flattened <- edges_flatten(adj_mat, ggplot.canvas$data$x, ggplot.canvas$data$y)
  ggplot.canvas + geom_segment(aes(x = x, y = y, xend = xend, yend = yend), size = 0.03, data = flattened)
}

plot_links <- function(ggplot.canvas, adj_mat){
  flattened <- edges_flatten(adj_mat, ggplot.canvas$data$x, ggplot.canvas$data$y)
  ggplot.canvas + geom_segment(aes(x = x, y = y, xend = xend, yend = yend), size = 0.03, arrow = arrow(length = unit(0.01, "npc")), data = flattened)
}

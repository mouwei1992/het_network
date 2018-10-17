# just to make it easy to navigate to my working folder
# not to be updated online
setwd("~/Dropbox/het_network_project/version6/src/experiments/real_data")

source('../../post_analysis/auc_extract.R')
source('../../post_analysis/edges_flatten.R')
source('../../post_analysis/network_plot.R')
source('../../post_analysis/pairwise_dist.R')

# reading adjacency matrices
cla_data_path <- "../../../data_file/classroom/klas12b-net-"
read_cla <- function(file_name) as.matrix(read.csv(
  paste(cla_data_path, file_name, sep = ""),sep = "", header = FALSE,
  blank.lines.skip = TRUE))

existing_nodes = c(1:20, 22:26)

cla_1_adj_t1 <- read_cla("1.dat")
cla_1_adj_t2 <- read_cla("2.dat")
cla_1_adj_t3 <- read_cla("3.dat")
cla_1_adj_t3 <- cla_1_adj_t3[existing_nodes, existing_nodes]
cla_1_adj_t4 <- read_cla("4.dat")
cla_1_adj_t4 <- cla_1_adj_t4[existing_nodes, existing_nodes]

# reading simulated positions
read_pos <- function(file_name) as.matrix(read.csv(
    paste("classroom_conf_1_position_result",file_name, sep=''), header = FALSE))


pos_1 <- read_pos("0.csv")
pos_2 <- read_pos("1.csv")
pos_3 <- read_pos("2.csv")
pos_3 <- pos_3[, existing_nodes]
pos_4 <- read_pos("3.csv")
pos_4 <- pos_4[, existing_nodes]

nodes.type <- rep("student", 26)
nodes.type2 <- rep(0, 25)

plot_labels <- function(ggplot.canvas){
  (
    ggplot.canvas + geom_text(aes(
      x = x,
      y = y,
      label = 1:26
    ))
  )
}

plot_labels_2 <- function(ggplot.canvas){
  (
    ggplot.canvas + geom_text(aes(
      x = x,
      y = y,
      label = existing_nodes)
  ))
}
# plotting..
classroom.canvas <- init_canvas(pos_1, nodes.type)
pts_canvas <- plot_points(classroom.canvas)
plot_labels(plot_links(pts_canvas, cla_1_adj_t1))

classroom.canvas <- init_canvas(pos_2, nodes.type)
pts_canvas <- plot_points(classroom.canvas)
plot_labels(plot_links(pts_canvas, cla_1_adj_t2))

classroom.canvas <- init_canvas(pos_3, nodes.type2)
pts_canvas <- plot_points(classroom.canvas)
plot_labels_2(plot_links(pts_canvas, cla_1_adj_t3))

classroom.canvas <- init_canvas(pos_4, nodes.type2)
pts_canvas <- plot_points(classroom.canvas)
plot_labels_2(plot_links(pts_canvas, cla_1_adj_t4))


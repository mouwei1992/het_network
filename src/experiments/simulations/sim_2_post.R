source_path = "~/Dropbox/het_network_project/version6/src/experiments/simulations"
setwd(source_path)

source("../../post_analysis/network_plot.R")
source("../../post_analysis/auc_extract.R")
source("../../post_analysis/edges_flatten.R")
# source("sim_1_readadj.R")
library(vegan)

read_pos <- function(file_name) as.matrix(read.csv(file_name, header = FALSE))
pos_1 <- read_pos("sim2_positions_sim1.csv")
pos_1_true <- read_pos("sim2_positions_true0.csv")

nodes.type = rep(c(rep(0, 8), rep(1, 25), rep(2, 47)), 3)

canv_sim_2_t1 <- init_canvas(pos_1, nodes.type)
plot_points(canv_sim_2_t1)

canv_sim_2_t1_true <- init_canvas(pos_1_true, nodes.type)
plot_points(canv_sim_2_t1_true)

rot <- procrustes(t(pos_1_true), t(pos_1), scale = F)$rotation
pos_1 <- rot %*% pos_1

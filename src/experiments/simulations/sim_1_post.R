source_path = "~/Dropbox/het_network_project/version6/src/experiments/simulations"
setwd(source_path)

source("../../post_analysis/network_plot.R")
source("../../post_analysis/auc_extract.R")
source("../../post_analysis/edges_flatten.R")
source("../../post_analysis/pairwise_dist.R")
source("sim_1_readadj.R")

read_pos <- function(file_name, sep = ',') as.matrix(read.csv(file_name, header = FALSE, sep = sep))

pos_1 <- read_pos("sim_1_position_result0.csv")
pos_2 <- read_pos("sim_1_position_result1.csv")
pos_3 <- read_pos("sim_1_position_result2.csv")
pos_4 <- read_pos("sim_1_position_result3.csv")

true_pos_1 <- read_pos("sim_1_true_position1.csv", ' ')
true_pos_2 <- read_pos("sim_1_true_position2.csv", ' ')
true_pos_3 <- read_pos("sim_1_true_position3.csv", ' ')
true_pos_4 <- read_pos("sim_1_true_position4.csv", ' ')

nodes.type <- c(rep(1,60), rep(2,90))

canv_sim_1_t1 <- init_canvas(pos_1, nodes.type)
pt_sim_1_t1 <- plot_points(canv_sim_1_t1)
# plot_edges(pt_sim_1_t1, sim_1_adj_t1)
pt_sim_1_t1
# plot_links(pt_sim_1_t1, sim_1_adj_t1)

canv_sim_1_t2 <- init_canvas(pos_2, nodes.type)
plot_points(canv_sim_1_t2)

canv_sim_1_t3 <- init_canvas(pos_3, nodes.type)
plot_points(canv_sim_1_t3)

canv_sim_1_t4 <- init_canvas(pos_4, nodes.type)
plot_points(canv_sim_1_t4)

dist_ratio_t1 = dist_ratio(pos_1, true_pos_1)
dist_ratio_t2 = dist_ratio(pos_2, true_pos_2)
dist_ratio_t3 = dist_ratio(pos_3, true_pos_3)
dist_ratio_t4 = dist_ratio(pos_4, true_pos_4)

dist_ratio = c(dist_ratio_t1, dist_ratio_t2, dist_ratio_t3, dist_ratio_t4)

qplot(dist_ratio, geom="histogram", binwidth = 0.01, xlim = c(0,5))

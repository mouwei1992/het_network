sim_1_data_path <- "../../../data_file/sim_1/sim_1_t"
read_sim1 <- function(file_name) read.csv(
  paste(sim_1_data_path, file_name, sep = ""), header = FALSE)

sim_1_adj_t1 <- read_sim1("1.csv")
sim_1_adj_t2 <- read_sim1("2.csv")
sim_1_adj_t3 <- read_sim1("3.csv")
sim_1_adj_t4 <- read_sim1("4.csv")
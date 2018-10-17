#include "../MCMC_module/dynamic_het.cpp"

int main(void){
  cluster_info c_info = read_c_info("sim2_pos_info.json");
  cluster_generator sim2_pos_generator(c_info);
  sim2_pos_generator.generate();

 network_structure n_str;
 network_model_pars n_pars;

 n_str.adj_matrix = Cube<Q_Q>(240, 240, 4, fill::zeros);
 n_str.node_type_each = Col<Q_Q>(240, fill::zeros);
 // 24 authors, 75 papers, 141 words
 n_str.node_type_each(span(8, 32)).fill(1);
 n_str.node_type_each(span(88, 112)).fill(1);
 n_str.node_type_each(span(168, 192)).fill(1);

 n_str.node_type_each(span(33, 79)).fill(2);
 n_str.node_type_each(span(113, 159)).fill(2);
 n_str.node_type_each(span(193, 239)).fill(2);

 // each paper only appear in one year

 n_str.regression_type = Mat<char>(3,3, fill::zeros);
 n_str.regression_type.fill(NOREG);
 n_str.regression_type(0,1) = LOGISTIC;
 n_str.regression_type(1,2) = POISSON;
 n_str.space_dim = 2;

 fill_structure_info(n_str);

 n_pars.positions = sim2_pos_generator.dup_positions();
 n_pars.betas = mat(3,3, fill::ones);
 n_pars.betas.at(0,1) = 0.8;
 n_pars.betas.at(1,2) = 0.2;
 n_pars.radius = mat(3,3, fill::zeros);
 n_pars.radius(0,1) = 0.3;
 n_pars.radius(1,2) = 0.7;
 dy_het_generator(n_str, n_pars, 0);

 network_info sim2_info;
 sim2_info.p_structure = n_str;

 sim2_info.p_prior = {
   1,
   1,
   0.5
 };

 sim2_info.p_update = {
   10000,
   0,
   2000,
   vec(3, fill::zeros),
   0.1,
   10000
 };

 sim2_info.p_update.positions_step_size.fill(0.2);
 dy_het_network sim2_net(sim2_info);

 sim2_net.update(10000);

 // true parameters
 save_result(n_pars, "", "../experiments/simulations/sim2_positions_true");
 save_result(sim2_net.duplicate_current_pars(), "", "../experiments/simulations/sim2_positions_sim");
 sim2_net.gen_log(n_pars, "../experiments/simulations/sim2_with_procrustes_per1000.log");

return 0;
}

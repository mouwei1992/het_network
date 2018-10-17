#include<string>

#include "../../utils/parse_config.h"
#include "../../utils/save_result.h"
#include "../../utils/save_fitted.h"
#include "../../utils/save_acceptance.h"
#include "../../utils/save_distance.h"
#include "../../MCMC_module/dynamic_het.cpp"

int main(){
  // set random seed
  arma_rng::set_seed_random();

  const string DATA_FILE_PATH = "../../../data_file/sim_1/";
  network_info sim1_info = parse_config("sim1_config.json", DATA_FILE_PATH);
  sim1_info.p_structure.regression_type = Mat<char>(2,2);
  sim1_info.p_structure.regression_type.fill(LOGISTIC);
  dy_het_network sim1_net(sim1_info);

  sim1_net.update(4000);

  // print parameters
  save_result(sim1_net.duplicate_current_pars());

  save_result(sim1_net.duplicate_current_pars(), "sim_1_par_result", "sim_1_position_result");
  save_distance(sim1_net.duplicate_current_cache().dist_matrix, "sim_1_dist_result");
  save_fitted(sim1_net.fitted_adj_mat(), "sim_1_fitted");
  save_acceptance(sim1_net, "sim1_accepted");
  return 0;
}

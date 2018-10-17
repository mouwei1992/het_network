#include <string>
#include "../MCMC_module/dynamic_het.cpp"
#include "../utils/parse_config.h"
#include "../utils/save_result.h"

int main(){
  string config_file_name = "../experiments/simulations/sim1_config.json";
  network_info n_info_sim1 = parse_config(config_file_name, "../../data_file/sim_1/");
  dy_het_network dhn_sim1(n_info_sim1);

  dhn_sim1.update(1000);
  print_pars(std::cout, dhn_sim1.duplicate_current_pars());
  save_result(dhn_sim1.duplicate_current_pars(), "sim_1_pars.csv", "sim_1_last_position");
  return 0;
}

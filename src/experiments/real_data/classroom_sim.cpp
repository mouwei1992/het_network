#include "../../utils/parse_config.h"
#include "../../utils/save_result.h"
#include "../../utils/save_fitted.h"
#include "../../utils/save_acceptance.h"
#include "../../MCMC_module/dynamic_het.cpp"

int main(){
  arma_rng::set_seed_random();
  const string DATA_FILE_PATH = "../../../data_file/classroom/";
  network_info classroom_info = parse_config("classroom_configure_1.json", DATA_FILE_PATH);
  dy_het_network classroom_net(classroom_info);

  classroom_net.update(10000);

  classroom_net.show_acceptance();
  save_result(classroom_net.duplicate_current_pars(), "classroom_conf_1_par_result", "classroom_conf_1_position_result");
  save_fitted(classroom_net.fitted_adj_mat(), "classroom_fitted");
  return 0;
}

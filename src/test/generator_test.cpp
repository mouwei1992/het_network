#include "../MCMC_module/dynamic_het.cpp"
#include "../utils/dy_het_generator.cpp"
#include "../utils/fill_info.h"

int main(void){
  network_structure n_str;
  network_model_pars n_pars;

  n_str.adj_matrix = Cube<Q_Q>(2, 2, 2, fill::zeros);
  n_str.node_type_each = Col<Q_Q>(2, fill::zeros);
  n_str.regression_type = Mat<char>(1,1, fill::zeros);
  n_str.space_dim = 2;
  n_str.regression_type(0,0) = LOGISTIC;


  fill_structure_info(n_str);


  n_pars.positions = cube(2, 2, 2, fill::randu);
  n_pars.betas = mat(1,1, fill::ones);
  n_pars.radius = mat(1,1, fill::ones);
  dy_het_generator(n_str, n_pars, 0);


  // testing procrustes_max_likelihood
  network_info test_info;
  test_info.p_structure = n_str;
  test_info.p_prior = {1, 1, 1};
  test_info.p_update = {1,1,1,mat(1,1,fill::ones),1,1};
  dy_het_network test_net(test_info);

  cout << test_net.duplicate_current_pars().positions;
  test_net.procrustes_max_likelihood();
  cout << test_net.duplicate_current_pars().positions;


  // it should work..
  return 0;
}

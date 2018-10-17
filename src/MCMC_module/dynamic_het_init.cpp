#ifndef DYNAMIC_HET_INIT
#define DYNAMIC_HET_INIT

dy_het_network::dy_het_network(network_info &n_info){

  fill_structure_info(n_info.p_structure);
  structure_pars = n_info.p_structure;
  update_pars = n_info.p_update;
  prior_pars = n_info.p_prior;

  init_draw_from_prior();
  init_buf();
  init_acceptance();
  draw_missing();
  update_current_dist();
  update_current_llh();
}

void dy_het_network::init_buf(void){
  const n_n &n_nodes = structure_pars.n_nodes;
  const n_n &t_length = structure_pars.time_length;

  current_buf.dist_matrix = cube(n_nodes,n_nodes,t_length);
  current_buf.log_likelihood_mat = cube(n_nodes,n_nodes,t_length);
  proposal_buf.dist_matrix = cube(n_nodes,n_nodes,t_length);
  proposal_buf.log_likelihood_mat = cube(n_nodes,n_nodes,t_length);

  current_buf.dist_matrix.zeros();
  proposal_buf.dist_matrix.zeros();
  current_buf.log_likelihood_mat.zeros();
  proposal_buf.log_likelihood_mat.zeros();
}

void dy_het_network::init_acceptance(void){
  const n_n &n_nodes = structure_pars.n_nodes;
  const n_n &t_length = structure_pars.time_length;
  const Q_Q &n_types = structure_pars.n_types;

  total_acceptance.positions_accepted = Mat<n_n>(n_nodes,t_length);
  total_acceptance.betas_accepted = Mat<n_n>(n_types,n_types);

  total_acceptance.positions_accepted.zeros();
  total_acceptance.betas_accepted.zeros();
  total_acceptance.radius_accepted = 0;
  total_acceptance.s2_accepted = 0;
  total_acceptance.t2_accepted = 0;
}

void dy_het_network::init_draw_from_prior(void){
  current_pars.s2 = r_hf_cauchy(prior_pars.s2_scale);
  current_pars.t2 = r_hf_cauchy(prior_pars.t2_scale);

  const Q_Q &n_types = structure_pars.n_types;

  current_pars.radius = mat(n_types,n_types);
  current_pars.radius.ones();
  current_pars.radius.for_each([this](double &elem) {elem *= update_pars.radius_kappa;});

  current_pars.betas = mat(n_types,n_types);
  for(Q_Q type_1 = 0; type_1 < n_types; type_1++){
    for(Q_Q type_2 = 0; type_2 < n_types; type_2++){
      current_pars.betas(type_1,type_2) = r_hf_cauchy(prior_pars.beta_scale);
      if(structure_pars.regression_type.at(type_1, type_2) == NOREG){
        current_pars.radius.at(type_1, type_2) = 0;
        current_pars.betas.at(type_1, type_2) = 0;
      }
    }
  }

  r_dirichlet(current_pars.radius);

  const n_n &t_length = structure_pars.time_length;
  const n_n &n_nodes = structure_pars.n_nodes;
  const Q_Q &universe_dimension = structure_pars.space_dim;
  current_pars.positions = cube(universe_dimension,n_nodes,t_length);
  current_pars.positions.slice(0).randn();
  current_pars.positions.slice(0).for_each([this](double &elem) {elem *= sqrt(current_pars.s2);});
  for(n_n t = 1; t < t_length; t++){
    current_pars.positions.slice(t).randn();
    current_pars.positions.slice(t).for_each([this](double &elem) {elem *= sqrt(current_pars.t2);});
    current_pars.positions.slice(t) += current_pars.positions.slice(t-1);
  }

  proposal_pars = current_pars;
  return;
}

#endif

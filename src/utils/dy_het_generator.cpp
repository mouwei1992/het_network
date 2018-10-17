
void dy_het_generator(network_structure &n_str, network_model_pars &n_pars, const double &p_missing){
  arma_rng::set_seed_random();

  // n_pars only need to provide positions info
  // n_str provides #nodes in each type
  // and regression type
  n_n &n_nodes = n_str.n_nodes;
  n_n &t_length = n_str.time_length;

  Cube<Q_Q> &adj_mat = n_str.adj_matrix;

  // inin adj mat
  cube &positions = n_pars.positions;
  t_length = positions.n_slices;
  n_nodes = positions.n_cols;
  adj_mat = Cube<Q_Q>(n_nodes, n_nodes, t_length);

  double ru = 0;
  double dist = 0;
  double eta;
  double glm_mean = 0;

  for(n_n node_1 = 0; node_1 < n_nodes; node_1++){
    Q_Q &node_1_type = n_str.node_type_each.at(node_1);
    for(n_n node_2 = 0; node_2 < n_nodes; node_2++){
      if(node_1 == node_2) continue;
      Q_Q &node_2_type = n_str.node_type_each.at(node_2);
      if(n_str.regression_type.at(node_1_type, node_2_type) == NOREG) continue;

      for(n_n t = 0; t < t_length; t++){
        ru = randu();
        if((ru < p_missing) && (n_str.regression_type.at(node_1_type, node_2_type) == LOGISTIC)){
          adj_mat.at(node_1, node_2, t) = MISSING;
          continue;
        }
        ru = randu();
        double &beta = n_pars.betas.at(node_1_type, node_2_type);
        double &radius = n_pars.radius.at(node_1_type, node_2_type);

        dist = norm(positions.slice(t).col(node_1) - positions.slice(t).col(node_2));
        eta = beta * (radius / dist - dist / radius);

        switch(n_str.regression_type.at(node_1_type, node_2_type)){
          case LOGISTIC:
            glm_mean = 1 / ( 1 + exp(-eta));
            if(ru < glm_mean)
              adj_mat.at(node_1, node_2, t) = 1;
            else adj_mat.at(node_1, node_2, t) = 0;
            break;
          case POISSON:
            glm_mean = eta>0? (eta+1):exp(eta);
            // glm_mean = eta > 0? eta:0;
            adj_mat.at(node_1, node_2, t) = r_poisson(glm_mean);
            break;
          default:
            cout << "pls specify link type" << endl;
            return;
        }
      }
    }
  }
  return;
}

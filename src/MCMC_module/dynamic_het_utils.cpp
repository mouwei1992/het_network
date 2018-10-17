#ifndef DYNAMIC_HET_UTILS
#define DYNAMIC_HET_UTILS

void dy_het_network::procrustes_max_likelihood(void){
  cube &positions = current_pars.positions;
  const n_n &n_nodes = structure_pars.n_nodes;
  const Mat<Q_Q> &rmed = structure_pars.node_removed;
  n_n t_length = structure_pars.time_length;

  uvec existing_nodes, common_nodes;
  vec last_center(structure_pars.space_dim, fill::zeros);
  vec current_center;
  mat rotation_mat;

  for(n_n t = 0; t < t_length; t++){
    existing_nodes = find(rmed.col(t) == 0);

    if(t > 0){
      common_nodes = find((rmed.col(t) == 0) && (rmed.col(t - 1) == 0));
      last_center = mean(positions.slice(t - 1).cols(common_nodes), 1);
      // find rotation mat
      rotation_mat = procrusted_rotate(
        positions.slice(t - 1).cols(common_nodes),
        positions.slice(t).cols(common_nodes)
      );
    }
    // when t is 0, just translate to origin
    else common_nodes = existing_nodes;

    current_center = mean(positions.slice(t).cols(common_nodes), 1);

    for(n_n node = 0; node < n_nodes; node++){
      if(rmed(node, t)) continue;
      positions.slice(t).col(node) -= current_center;
      if(t > 0) positions.slice(t).col(node) = rotation_mat * positions.slice(t).col(node);
      positions.slice(t).col(node) += last_center;
    }
  }
  return;
}

void dy_het_network::update_current_dist(void){
  cube &positions = current_pars.positions;
  cube &dist_mat = current_buf.dist_matrix;
  const n_n &n_nodes = structure_pars.n_nodes;
  const Cube<Q_Q> &adj_mat = structure_pars.adj_matrix;
  Mat<Q_Q> &rmed = structure_pars.node_removed;

  for(Q_Q t = 0; t < structure_pars.time_length; t++){
    for(n_n node_1 = 0; node_1 < n_nodes; node_1++){
      for(n_n node_2 = node_1; node_2 < n_nodes; node_2++){
        if(node_1 == node_2){
          // update distance between node_1,t and node_1,t-1
          if(rmed.at(node_1,t)) continue;
          if(t == 0) dist_mat(node_1,node_1,t) = norm(positions.slice(t).col(node_1));
          // though we generated some random positions
          // only use them when node's not removed
          else if(!rmed.at(node_1,t-1)) dist_mat(node_1,node_1,t) = norm(positions.slice(t).col(node_1) - positions.slice(t-1).col(node_1));
          else continue;
          dist_mat(node_1,node_1,t) *= dist_mat(node_1,node_1,t);
        }else{
          if(rmed.at(node_1, t) || rmed.at(node_2,t)) continue;
          // if(adj_mat.at(node_1,node_2,t) == REMOVED) continue;
          // node_1 < node_2
          dist_mat(node_1,node_2,t) = norm(positions.slice(t).col(node_1) - positions.slice(t).col(node_2));
        }
      }
    }
  }

  return;
}

inline double dy_het_network::llh_cal(Q_Q y, double &eta, const reg_type &r_type){
  double llh;
  switch(r_type){
    case LOGISTIC:
      if(eta > 25) llh = -eta;
      else if(eta < -25) llh = 0;
      else llh = -log(exp(float(eta)) + 1);
      // if(y) llh += eta;
      switch(y){
        case UNLINKED:
        case MISSING_0:
          break;
        case LINKED:
        case MISSING_1:
          llh += eta;
          break;
        default:
          cout << "unrecognized link type of "<< int(y) << endl;
          break;
      }
      break;
    case POISSON:
      if(eta < 0 ) llh = - exp(eta) + y * eta - lgamma(y + 1);
      else llh = - (eta + 1) + y * log(eta + 1) - lgamma(y + 1);
      break;
    case TEST:
      llh = eta;
      break;
    // a faster but less accurate version
    case LOGISTIC_F:
      llh = logistic_prox(y,eta);
      break;

    case NOREG:
      cout << "should be calculating this" << endl;
      break;

    default:
      cout  << "unrecognized link type of "<< int(y) << endl;
      break;
  }
  return llh;
}

void *dy_het_network::update_current_llh_single_helper(void *context){
    thread_info *th_if = (thread_info *)context;
    ( ((dy_het_network *)(th_if->network_pt)) ->update_current_llh_single(th_if->thread_id));
    return 0;
}

void dy_het_network::update_current_llh(void){
  pthread_t tid1, tid2, tid3, tid4;
  thread_info t1_ = {this,0};
  thread_info t2_ = {this,1};
  thread_info t3_ = {this,2};
  thread_info t4_ = {this,3};
  pthread_create(&tid1, NULL, &dy_het_network::update_current_llh_single_helper, &t1_);
  pthread_create(&tid2, NULL, &dy_het_network::update_current_llh_single_helper, &t2_);
  pthread_create(&tid3, NULL, &dy_het_network::update_current_llh_single_helper, &t3_);
  pthread_create(&tid4, NULL, &dy_het_network::update_current_llh_single_helper, &t4_);
  pthread_join(tid1,NULL);
  pthread_join(tid2,NULL);
  pthread_join(tid3,NULL);
  pthread_join(tid4,NULL);
  return;
}

void dy_het_network::update_current_llh_single(Q_Q &n_th){
  cube &dist_mat = current_buf.dist_matrix;
  cube &llh_mat = current_buf.log_likelihood_mat;
  n_n &n_nodes = structure_pars.n_nodes;
  const Mat<Q_Q> &rmed = structure_pars.node_removed;
  Mat<char> r_type = structure_pars.regression_type;
  Q_Q node_1_type;
  Q_Q node_2_type;
  double dist;
  double eta;
  for(n_n node_1 = 0; node_1 < n_nodes; node_1++){
    if(node_1 % 4 != n_th) continue;
    node_1_type = structure_pars.node_type_each[node_1];
    for(n_n node_2 = 0; node_2 < n_nodes; node_2++){
      if(node_1 == node_2) continue;
      node_2_type = structure_pars.node_type_each[node_2];
      // if it is not what we are interested
      // if(!is_valid_link(node_1_type, node_2_type)) continue;
      if(r_type.at(node_1_type, node_2_type) == NOREG) continue;
      for(Q_Q t = 0; t < structure_pars.time_length; t++){
        if(rmed.at(node_1, t) || rmed.at(node_2,t)) continue;
        // if(structure_pars.adj_matrix.at(node_1,node_2,t) == REMOVED) continue;
        const double &dist = node_1 < node_2 ? dist_mat.at(node_1,node_2,t):dist_mat.at(node_2,node_1,t);
        eta = eta_cal(dist, current_pars.betas.at(node_1_type,node_2_type), current_pars.radius.at(node_1_type,node_2_type));
        llh_mat.at(node_1,node_2,t) = llh_cal(structure_pars.adj_matrix.at(node_1,node_2,t), eta, (reg_type)r_type.at(node_1_type, node_2_type));
      }
    }
  }
  return;
}

void dy_het_network::draw_missing(void){
  // currently single threaded
  double eta = 0;
  double ru = 0;
  n_n &n_nodes = structure_pars.n_nodes;
  n_n &t_length = structure_pars.time_length;
  Cube<Q_Q> &adj_mat = structure_pars.adj_matrix;
  const Mat<Q_Q> &rmed = structure_pars.node_removed;
  cube &current_dist = current_buf.dist_matrix;
  Mat<char> r_type = structure_pars.regression_type;

  for(n_n node_1 = 0; node_1 < n_nodes; node_1++){
    for(n_n t = 0; t < t_length; t++){
      if(rmed.at(node_1,t)) continue;
      Q_Q &node_1_type = structure_pars.node_type_each(node_1);
      for(n_n node_2 = 0; node_2 < n_nodes; node_2++){
        // check if the edge is missing
        // missing: 9/ - 91/ - 90
        // deleted 10
        if((adj_mat.at(node_1, node_2, t) < 0)||(adj_mat.at(node_1, node_2, t) == MISSING)){
          Q_Q &node_2_type = structure_pars.node_type_each(node_2);
          if(r_type.at(node_1_type, node_2_type) != LOGISTIC) continue;
          // if(!is_valid_link(node_1_type, node_2_type)) continue;
          const double &dist = node_1 < node_2 ?
            current_dist.at(node_1,node_2,t):current_dist.at(node_2,node_1,t);
          ru = randu();
          eta = eta_cal(
            dist,
            current_pars.betas.at(node_1_type,node_2_type),
            current_pars.radius.at(node_1_type,node_2_type)
          );
          if( log(ru) < llh_cal(0, eta, (reg_type)r_type.at(node_1_type, node_2_type)) ) // with prob pr(unlinked | eta)
            adj_mat.at(node_1,node_2,t) = MISSING_0;
          else adj_mat.at(node_1,node_2,t) = MISSING_1;
        }
      }
    }
  }
  return;
}

// default print on console
void dy_het_network::show_acceptance(ostream& ost) const{
  ost << "positions accepted:" << endl << total_acceptance.positions_accepted;
  ost << "beta accepted:" << endl << total_acceptance.betas_accepted;
  ost << "radius accepted: " << total_acceptance.radius_accepted << endl;
  return;
}

void dy_het_network::save_acceptance(const string &file_start = "") const{
  if(file_start.empty()) show_acceptance(cout);
  else{
    ofstream file_writer;
    // warning: will overwrite data
    file_writer.open(file_start, std::ofstream::out | std::ofstream::trunc);
    show_acceptance(file_writer);
  }
}

cube dy_het_network::fitted_adj_mat(void) const{
  const n_n &n_nodes = structure_pars.n_nodes;
  const n_n &t_length = structure_pars.time_length;
  cube fitted_adj_matrix(n_nodes, n_nodes, t_length, fill::zeros);
  const Mat<Q_Q> &rmed = structure_pars.node_removed;
  const cube &dist_mat = current_buf.dist_matrix;

  double eta;

  for(n_n node_1 = 0; node_1 < n_nodes; node_1++){
    const Q_Q &node_1_type = structure_pars.node_type_each(node_1);
    for(n_n node_2 = 0; node_2 < n_nodes; node_2++){
      if(node_1 == node_2) continue;
      const Q_Q &node_2_type = structure_pars.node_type_each(node_2);
      // if(!is_valid_link(node_1_type, node_2_type)) continue;
      for(n_n t = 0; t < t_length; t++){
        if(rmed.at(node_1, t) || rmed.at(node_2,t)) continue;
        // if(structure_pars.adj_matrix.at(node_1, node_2, t) == REMOVED) continue;
        const double &dist = node_1 < node_2 ? dist_mat.at(node_1,node_2,t):dist_mat.at(node_2,node_1,t);
        const double &beta = current_pars.betas.at(node_1_type, node_2_type);
        const double &radius = current_pars.radius.at(node_1_type, node_2_type);
        eta = eta_cal(dist, beta, radius);
        // prob of link
        switch(structure_pars.regression_type.at(node_1_type, node_2_type)){
          case LOGISTIC:
            fitted_adj_matrix.at(node_1, node_2, t) = 1/(1 + exp(-eta));
            break;

          case POISSON:
            fitted_adj_matrix.at(node_1, node_2, t) = exp(eta);
            break;

          case NOREG:
            break;

          default:
            cout << "pls specify regression type" << endl;
            break;
        }
      }
    }
  }
  return fitted_adj_matrix;
}

network_model_pars dy_het_network::duplicate_current_pars(void) const{
  return current_pars;
}

network_cache dy_het_network::duplicate_current_cache(void) const{
  return current_buf;
}

void dy_het_network::gen_log(const network_model_pars &true_pars, const string &log_file_name = ""){
  // save true parameters, required
  save_result(true_pars, log_file_name);
  save_result(current_pars, log_file_name);
  // save prior parameters
  save_prior(prior_pars, log_file_name);
  // save MCMC setting
  save_stepsize(update_pars, log_file_name);

  return;
}

#endif

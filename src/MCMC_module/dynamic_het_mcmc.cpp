#ifndef DYNAMIC_HET_MCMC
#define DYNAMIC_HET_MCMC

int dy_het_network::update(n_n steps){

  cout << "Burn-in steps: " << update_pars.burn_in_steps << endl;
    // main MCMC update step

  for(int step = 0; step < steps; step++){
    update_pars.current_step++;

    if(update_pars.current_step == update_pars.burn_in_steps)
      cout << "Burn-in complets." << endl;

    // won't allow for predefined total steps, but it actually is redundant
    if(update_pars.current_step == update_pars.total_steps){
      cout << "Sample complete." << endl;
      return 0;
    }

    if(update_pars.current_step % 100 == 0)
      cout << "Current Step: " << update_pars.current_step << "/" << update_pars.total_steps << endl;

    draw_missing();
    update_positions();
    update_s2();
    update_t2();
    update_betas();
    update_radius();
    if((step > update_pars.burn_in_steps) && (step % 1000 == 0)) procrustes_max_likelihood();
  }
  return 1;
}

void dy_het_network::update_positions(void){
  const n_n &n_nodes = structure_pars.n_nodes;
  const n_n &t_length = structure_pars.time_length;
  // update nodes in a random sequence
  Col<n_n> update_seq = shuffle(linspace<Col<n_n>> (0,n_nodes - 1,n_nodes));
  for(n_n node_cnt = 0; node_cnt < n_nodes; node_cnt++){
    for(n_n t = 0; t < t_length; t++) update_position(update_seq(node_cnt),t);
  }
  return;
}

void dy_het_network::update_position(const n_n &node, const n_n &t){

  const Mat<Q_Q> &rmed = structure_pars.node_removed;
  if(rmed.at(node,t)) return;

  const n_n &n_nodes = structure_pars.n_nodes;
  const Q_Q &node_type = structure_pars.node_type_each(node);
  mat &current_nodes_positions = current_pars.positions.slice(t);
  // actually i am not sure what i am doing..
  vec current_position = current_nodes_positions.unsafe_col(node);
  // same here
  vec proposal_position = proposal_pars.positions.slice(t).unsafe_col(node);
  Mat<char> r_type = structure_pars.regression_type;
  cube &current_dist = current_buf.dist_matrix;
  cube &current_llh = current_buf.log_likelihood_mat;
  cube &proposal_dist = proposal_buf.dist_matrix;
  cube &proposal_llh = proposal_buf.log_likelihood_mat;
  double dist;
  double eta;
  double accept_prob = 0;
  double ru = randu();

  proposal_position.randn();
  proposal_position.for_each([this, &node_type](double &elem) {elem *= (update_pars.positions_step_size)[node_type];});
  proposal_position += current_position;

  for(n_n node_2 = 0; node_2 < n_nodes; node_2++){
    if(rmed.at(node_2, t)) continue;
    // if(structure_pars.adj_matrix.at(node,node_2,t) == REMOVED) continue;
    if(node == node_2) continue;
    Q_Q &node_2_type = structure_pars.node_type_each[node_2];

    // no need for distances of no interest

    accept_prob -= (current_llh.at(node,node_2,t) + current_llh.at(node_2,node,t));
    dist = norm(proposal_position - current_nodes_positions.col(node_2));
    if(node < node_2) proposal_dist.at(node,node_2,t) = dist;
    else proposal_dist.at(node_2,node,t) = dist;

    // node -> node_2
    if(r_type.at(node_type, node_2_type) != NOREG){
      eta = eta_cal(dist,current_pars.betas.at(node_type,node_2_type),current_pars.radius.at(node_type,node_2_type));
      proposal_llh.at(node,node_2,t) = llh_cal(structure_pars.adj_matrix.at(node,node_2,t),eta, (reg_type)r_type.at(node_type, node_2_type));
    }

    // the other way
    if(r_type.at(node_2_type, node_type) != NOREG){
      eta = eta_cal(dist,current_pars.betas.at(node_2_type,node_type),current_pars.radius.at(node_2_type,node_type));
      proposal_llh.at(node_2,node,t) = llh_cal(structure_pars.adj_matrix.at(node_2,node,t),eta, (reg_type)r_type.at(node_2_type, node_type));
    }

    accept_prob += (proposal_llh.at(node,node_2,t) + proposal_llh.at(node_2,node,t));
  }

  if(t == 0){
    dist = dot(proposal_position,proposal_position);
    proposal_dist.at(node,node,0) = dist;
    accept_prob -= dist / (2 * current_pars.s2);
    double &old_dist = current_dist(node,node,0);
    accept_prob += old_dist / (2 * current_pars.s2);
  }

  // if the node exists at t-1
  else if(!rmed.at(node, t - 1)){
    vec last_position = current_pars.positions.slice(t - 1).unsafe_col(node);
    // it is redundant to do it here but it is not in the main loop..
    vec diff_position = proposal_position - last_position;
    dist = dot(diff_position, diff_position);
    proposal_dist.at(node,node,t) = dist;
    accept_prob -= dist * dist / (2 * current_pars.t2);
    double &old_dist = current_dist.at(node,node,t);
    accept_prob += old_dist * old_dist / (2 * current_pars.t2);
  }

  // if t is not T, and the node is one that exists at (t + 1)
  if((t < structure_pars.time_length - 1)&&(!rmed.at(node, t+1))){
    // dangerous!
    vec next_position = current_pars.positions.slice(t + 1).unsafe_col(node);
    vec diff_position = proposal_position - next_position;
    dist = dot(diff_position, diff_position);
    proposal_dist.at(node,node,t + 1) = dist;
    accept_prob -= dist / (2 * current_pars.t2);
    double &old_dist = current_dist.at(node,node,t + 1);
    accept_prob += old_dist /(2 * current_pars.t2);
  }
  accept_prob = accept_prob > 0 ? 1:exp(accept_prob);

  if(ru < accept_prob){
    current_position = proposal_position;
    // update the buffers..
    current_dist.at(node,node,t) = proposal_dist.at(node,node,t);
    if((t < structure_pars.time_length - 1)&&(!rmed.at(node, t+1)))
     current_dist(node,node,t + 1) = proposal_dist.at(node,node,t + 1);
    for(n_n node_2 = 0; node_2 < n_nodes; node_2++){
      if(node == node_2) continue;
      if(rmed.at(node_2,t)) continue;
      // if(structure_pars.adj_matrix.at(node,node_2,t) == REMOVED) continue;
      current_llh.at(node,node_2,t) = proposal_llh.at(node,node_2,t);
      current_llh.at(node_2,node,t) = proposal_llh.at(node_2,node,t);
      if(node < node_2) current_dist.at(node,node_2,t) = proposal_dist.at(node,node_2,t);
      else current_dist.at(node_2,node,t) = proposal_dist.at(node_2,node,t);
    }
    if(update_pars.current_step >= update_pars.burn_in_steps) total_acceptance.positions_accepted.at(node,t)++;
  }
  return;
}

void dy_het_network::update_t2(void){
  double &new_t2 = proposal_pars.t2;
  double &old_t2 = current_pars.t2;
  n_n &aval_nodes = structure_pars.available_nodes_sum;
  double ru = randu();
  double accept_prob = 0;
  const n_n &n_nodes = structure_pars.n_nodes;
  const n_n &t_length = structure_pars.time_length;
  new_t2 = r_hf_cauchy(prior_pars.t2_scale);
  double step_sum = 0;
  for(n_n t = 1; t < t_length; t++)
    step_sum += trace(current_buf.dist_matrix.slice(t));
  accept_prob = aval_nodes * log(old_t2/new_t2) + step_sum/old_t2 - step_sum/new_t2;
  accept_prob /= 2.;
  accept_prob = accept_prob > 0 ? 1:exp(accept_prob);
  if(ru < accept_prob){
    old_t2 = new_t2;
    total_acceptance.t2_accepted++;
  }
  return;
}

void dy_het_network::update_s2(void){
  const n_n &n_nodes = structure_pars.n_nodes;
  double &new_s2 = proposal_pars.s2;
  double &old_s2 = current_pars.s2;
  double ru = randu();
  double accept_prob = 0;
  new_s2 = r_hf_cauchy(prior_pars.s2_scale);
  double step_sum = trace(current_buf.dist_matrix.slice(0));
  accept_prob = n_nodes * log(old_s2 / new_s2) + step_sum/old_s2 - step_sum/new_s2;
  accept_prob /= 2.;
  accept_prob = accept_prob > 0 ? 1:exp(accept_prob);
  if(ru < accept_prob){
    old_s2 = new_s2;
    total_acceptance.s2_accepted++;
  }
  return;
}

void dy_het_network::update_betas(void){
  const Q_Q &n_types = structure_pars.n_types;
  Mat<char> r_type = structure_pars.regression_type;

  for(Q_Q type_1 = 0; type_1 < n_types; type_1++){
    for(Q_Q type_2 = 0; type_2 < n_types; type_2++)
      if(r_type.at(type_1, type_2) != NOREG)  update_beta(type_1,type_2);
  }
  return;
}

void *dy_het_network::update_beta_single_helper(void *context){
  thread_info_b *th_if_b = (thread_info_b *) context;
  ((dy_het_network *)(th_if_b->network_pt))->update_beta_single(th_if_b->type_1, th_if_b -> type_2, th_if_b -> thread_id);
  return 0;
}

void dy_het_network::update_beta_single(const Q_Q &type_1, const Q_Q &type_2, const Q_Q &n_th){
  // only used to update the new llh_mat
  const cube &dist_mat = current_buf.dist_matrix;
  cube &new_llh_mat = proposal_buf.log_likelihood_mat;
  const n_n t_length = structure_pars.time_length;
  double &new_beta = proposal_pars.betas.at(type_1,type_2);
  const double &corr_radius = current_pars.radius.at(type_1,type_2);
  const Mat<Q_Q> &rmed = structure_pars.node_removed;
  const reg_type r_type = (reg_type)structure_pars.regression_type.at(type_1, type_2);

  double eta;

  const n_n &type_1_begin = structure_pars.each_type_start[type_1];
  const n_n &type_1_end = structure_pars.each_type_start[type_1 + 1];
  const n_n &type_2_begin = structure_pars.each_type_start[type_2];
  const n_n &type_2_end = structure_pars.each_type_start[type_2 + 1];

  for(n_n type_1_node = type_1_begin; type_1_node < type_1_end; type_1_node++){
    n_n &node_1 = structure_pars.nodes_rearrange[type_1_node];
    if(node_1 % 4 != n_th) continue;
    for(n_n type_2_node = type_2_begin; type_2_node < type_2_end; type_2_node++){
      n_n &node_2 = structure_pars.nodes_rearrange[type_2_node];
      if(node_1 == node_2) continue;
      for(n_n t = 0; t < t_length; t++){
        if(rmed.at(node_1, t) || rmed.at(node_2,t)) continue;
        // if(structure_pars.adj_matrix.at(node_1,node_2,t) == REMOVED) continue;
        const double &dist = node_1 < node_2 ? dist_mat.at(node_1,node_2,t):dist_mat.at(node_2,node_1,t);
        double &new_llh = new_llh_mat.at(node_1,node_2,t);
        eta = eta_cal(dist,new_beta,corr_radius);
        new_llh = llh_cal(structure_pars.adj_matrix.at(node_1,node_2,t), eta, r_type);
      }
    }
  }
  return;
}

void dy_het_network::update_beta(const Q_Q &type_1, const Q_Q &type_2){
  cube &old_llh_mat = current_buf.log_likelihood_mat;
  cube &new_llh_mat = proposal_buf.log_likelihood_mat;
  const n_n t_length = structure_pars.time_length;
  double &old_beta = current_pars.betas.at(type_1,type_2);
  double &new_beta = proposal_pars.betas.at(type_1,type_2);
  double ru = randu();
  const Mat<Q_Q> &rmed = structure_pars.node_removed;

  double accept_prob = 0;

  new_beta = (2 * ru - 1) * update_pars.beta_step_size + 1;
  if( new_beta * (1 + update_pars.beta_step_size) < 1) return;
  // transition prob..
  accept_prob = - log(new_beta);
  new_beta *= old_beta;

  pthread_t tid1, tid2, tid3, tid4;
  thread_info_b t1_ = {this,0, type_1, type_2};
  thread_info_b t2_ = {this,1, type_1, type_2};
  thread_info_b t3_ = {this,2, type_1, type_2};
  thread_info_b t4_ = {this,3, type_1, type_2};

  pthread_create(&tid1, NULL, &dy_het_network::update_beta_single_helper, &t1_);
  pthread_create(&tid2, NULL, &dy_het_network::update_beta_single_helper, &t2_);
  pthread_create(&tid3, NULL, &dy_het_network::update_beta_single_helper, &t3_);
  pthread_create(&tid4, NULL, &dy_het_network::update_beta_single_helper, &t4_);

  pthread_join(tid1,NULL);
  pthread_join(tid2,NULL);
  pthread_join(tid3,NULL);
  pthread_join(tid4,NULL);

  const n_n &type_1_begin = structure_pars.each_type_start[type_1];
  const n_n &type_1_end = structure_pars.each_type_start[type_1 + 1];
  const n_n &type_2_begin = structure_pars.each_type_start[type_2];
  const n_n &type_2_end = structure_pars.each_type_start[type_2 + 1];

  for(n_n type_1_node = type_1_begin; type_1_node < type_1_end; type_1_node++){
    n_n &node_1 = structure_pars.nodes_rearrange[type_1_node];
    for(n_n type_2_node = type_2_begin; type_2_node < type_2_end; type_2_node++){
      n_n &node_2 = structure_pars.nodes_rearrange[type_2_node];
      if(node_1 == node_2) continue;
      for(n_n t = 0; t < t_length; t++){
        if(rmed.at(node_1, t) || rmed.at(node_2,t)) continue;
        // if(structure_pars.adj_matrix.at(node_1,node_2,t) == REMOVED) continue;

        double &old_llh = old_llh_mat.at(node_1,node_2,t);
        double &new_llh = new_llh_mat.at(node_1,node_2,t);

        accept_prob += new_llh;
        accept_prob -= old_llh;
      }
    }
  }

  // times prior
  double &b_prior = prior_pars.beta_scale;
  accept_prob += log(1 + old_beta * old_beta / (b_prior * b_prior));
  accept_prob -= log(1 + new_beta * new_beta / (b_prior * b_prior));

  accept_prob = accept_prob < 0? exp(accept_prob):1;
  ru = randu();
  if(ru < accept_prob){
    old_beta = new_beta;
    for(n_n type_1_node = type_1_begin; type_1_node < type_1_end; type_1_node++){
      n_n &node_1 = structure_pars.nodes_rearrange[type_1_node];
      for(n_n type_2_node = type_2_begin; type_2_node < type_2_end; type_2_node++){
        if(type_1_node == type_2_node) continue;
        n_n &node_2 = structure_pars.nodes_rearrange[type_2_node];
        old_llh_mat.tube(node_1,node_2) = new_llh_mat.tube(node_1,node_2);
      }
    }
    if(update_pars.current_step >= update_pars.burn_in_steps) total_acceptance.betas_accepted.at(type_1,type_2)++;
  }
  return;
}


void dy_het_network::update_radius(void){
  double &kappa = update_pars.radius_kappa;
  mat &radius = current_pars.radius;
  mat old_radius = current_pars.radius;
  mat &new_radius = proposal_pars.radius;
  n_n accept = 0;
  double accept_prob = 0;
  double rand_u = randu();
  cube old_llh = current_buf.log_likelihood_mat;

  new_radius = old_radius;
  new_radius.for_each([&kappa](double &elem) {elem *= kappa;});
  r_dirichlet(new_radius);

  accept_prob -= accu(old_llh);
  radius = new_radius;
  update_current_llh();
  accept_prob += accu(current_buf.log_likelihood_mat);

  accept_prob = accept_prob>0 ? 1:exp(accept_prob);
  if(rand_u < accept_prob){
    total_acceptance.radius_accepted += 1;
  }else{
    radius = old_radius;
  }
  return;
}

#endif

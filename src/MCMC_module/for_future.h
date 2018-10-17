
/*
void dy_het_network::draw_missing(void){
  double eta = 0;
  double ru = 0;
  Cube<Q_Q> &adj_mat = structure_pars.adj_matrix;
  for(n_n node = 0; node < structure_pars.n_nodes; node++){
    for(n_n t = 0; t < structure_pars.time_length; t++){
      if(missing.at(node,t)){
        for(n_n node_2 = 0; node_2 < structure_pars.n_nodes; node_2++){
          if(missing.at(node_2,t) && (node_2 < node)) continue;
          const double &dist = node_1 < node_2 ? dist_mat.at(node_1,node_2,t):dist_mat.at(node_2,node_1,t);
          eta = eta_cal(dist, current_pars.betas.at(node_1_type,node_2_type), current_pars.radius.at(node_1_type,node_2_type));
          ru = randu();
          adj_mat.at(node,node_2,t) = ( log(ru) < llh_cal(0, eta) );
          eta = eta_cal(dist, current_pars.betas.at(node_2_type,node_1_type), current_pars.radius.at(node_2_type,node_1_type));
          ru = randu();
          adj_mat.at(node_2,node,t) = ( log(ru) < llh_cal(0, eta) );
        }
      }
    }
  }

  return;
}

*/

/*
void dy_het_network::update_positions(void){
  const n_n &n_nodes = structure_pars.n_nodes;
  const n_n &t_length = structure_pars.time_length;
  Col<n_n> update_seq = shuffle(linspace<Col<n_n>> (0,n_nodes - 1,n_nodes));

  cube &current_dist = current_buf.dist_matrix;
  cube &current_llh = current_buf.log_likelihood_mat;
  cube &proposal_dist = proposal_buf.dist_matrix;
  cube &proposal_llh = proposal_buf.log_likelihood_mat;

  double dist;
  double ru;
  double accept_prob;

  pthread_t tid1, tid2, tid3, tid4;
  pthread_mutex_t lk;
  pthread_cond_t position_proposed, llh_updated;

  pthread_mutex_init(&lk, NULL);
  pthread_cond_init(&position_proposed,NULL);
  pthread_cond_init(&llh_updated, NULL);

  thread_info_p t1_ = {this,0, &update_seq, &position_proposed, &llh_updated, &lk};
  thread_info_p t2_ = {this,1, &update_seq, &position_proposed, &llh_updated, &lk};
  thread_info_p t3_ = {this,2, &update_seq, &position_proposed, &llh_updated, &lk};
  thread_info_p t4_ = {this,3, &update_seq, &position_proposed, &llh_updated, &lk};

  pthread_create(&tid1, NULL, &dy_het_network::update_position_single_helper, &t1_);
  pthread_create(&tid2, NULL, &dy_het_network::update_position_single_helper, &t2_);
  pthread_create(&tid3, NULL, &dy_het_network::update_position_single_helper, &t3_);
  pthread_create(&tid4, NULL, &dy_het_network::update_position_single_helper, &t4_);

  for(n_n node_cnt = 0; node_cnt < n_nodes; node_cnt++){
    n_n &node = update_seq(node_cnt);
    for(n_n t = 0; t < t_length; t++){
      accept_prob = 0;
      ru = randu();
      Q_Q &node_type = structure_pars.node_type_each(node);
      mat &current_nodes_positions = current_pars.positions.slice(t);
      // actually i am not sure what i am doing..
      vec current_position = current_nodes_positions.unsafe_col(node);
      // same here
      vec proposal_position = proposal_pars.positions.slice(t).unsafe_col(node);

      proposal_position.randn();
      proposal_position.for_each([this, &node_type](double &elem) {elem *= (update_pars.positions_step_size)[node_type];});
      proposal_position += current_position;

      pthread_mutex_lock(&lk);
      pthread_mutex_lock(&lk);
      pthread_mutex_lock(&lk);
      pthread_mutex_lock(&lk);
      pthread_cond_broadcast(&position_proposed);

      if(t == 0){
        dist = dot(proposal_position,proposal_position);
        proposal_dist.at(node,node,0) = dist;
        accept_prob -= dist / (2 * current_pars.s2);
        double &old_dist = current_dist(node,node,0);
        accept_prob += old_dist / (2 * current_pars.s2);
      }else{
        // dangerous..
        vec last_position = current_pars.positions.slice(t - 1).unsafe_col(node);
        // it is redundant to do it here but it is not in the main loop..
        vec diff_position = proposal_position - last_position;
        dist = dot(diff_position, diff_position);
        proposal_dist.at(node,node,t) = dist;
        accept_prob -= dist * dist / (2 * current_pars.t2);
        double &old_dist = current_dist.at(node,node,t);
        accept_prob += old_dist * old_dist / (2 * current_pars.t2);
      }

      if(t < structure_pars.time_length - 1){
        // dangerous!
        vec next_position = current_pars.positions.slice(t + 1).unsafe_col(node);
        vec diff_position = proposal_position - next_position;
        dist = dot(diff_position, diff_position);
        proposal_dist.at(node,node,t + 1) = dist;
        accept_prob -= dist / (2 * current_pars.t2);
        double &old_dist = current_dist.at(node,node,t + 1);
        accept_prob += old_dist /(2 * current_pars.t2);
      }

      pthread_cond_wait(&llh_updated, &lk);
      pthread_cond_wait(&llh_updated, &lk);
      pthread_cond_wait(&llh_updated, &lk);
      pthread_cond_wait(&llh_updated, &lk);

      for(n_n node_2 = 0; node_2 < n_nodes; node_2++){
        if(node == node_2) continue;
        accept_prob -= (current_llh.at(node,node_2,t) + current_llh.at(node_2,node,t));
        accept_prob += (proposal_llh.at(node,node_2,t) + proposal_llh.at(node_2,node,t));
      }

      accept_prob = accept_prob > 0 ? 1:exp(accept_prob);

      if(ru < accept_prob){
        current_position = proposal_position;
        // update the buffers..
        current_dist.at(node,node,t) = proposal_dist.at(node,node,t);
        if(t < structure_pars.time_length - 1) current_dist(node,node,t + 1) = proposal_dist.at(node,node,t + 1);
        for(n_n node_2 = 0; node_2 < n_nodes; node_2++){
          current_llh.at(node,node_2,t) = proposal_llh.at(node,node_2,t);
          current_llh.at(node_2,node,t) = proposal_llh.at(node_2,node,t);
          if(node < node_2) current_dist.at(node,node_2,t) = proposal_dist.at(node,node_2,t);
          else current_dist.at(node_2,node,t) = proposal_dist.at(node_2,node,t);
        }
        if(update_pars.current_step >= update_pars.burn_in_steps) total_acceptance.positions_accepted.at(node,t) ++;
      }
    }
  }


  pthread_join(tid1,NULL);
  pthread_join(tid2,NULL);
  pthread_join(tid3,NULL);
  pthread_join(tid4,NULL);

  pthread_cond_destroy(&llh_updated);

  pthread_mutex_destroy(&lk);
  return;
}

void *dy_het_network::update_position_single_helper(void *context){
  thread_info_p *th_if_p = (thread_info_p *)context;
  ((dy_het_network *)(th_if_p -> network_pt)) -> update_position_single(th_if_p->update_seq , th_if_p -> thread_id, th_if_p -> position_proposed, th_if_p -> llh_updated, th_if_p -> lock);
  return 0;
}

void dy_het_network::update_position_single(Col<n_n> *update_seq, Q_Q &th_id, pthread_cond_t *position_proposed, pthread_cond_t *llh_updated, pthread_mutex_t *lock){
  const n_n &n_nodes = structure_pars.n_nodes;
  const n_n &t_length = structure_pars.time_length;

  cube &proposal_dist = proposal_buf.dist_matrix;
  cube &proposal_llh = proposal_buf.log_likelihood_mat;


  double dist;
  double eta;

  for(n_n node_cnt = 0; node_cnt < n_nodes; node_cnt++){
    n_n &node = (*update_seq)[node_cnt];
    Q_Q &node_type = structure_pars.node_type_each[node];

    for(n_n t = 0; t < t_length; t++){
      vec proposal_position = proposal_pars.positions.slice(t).unsafe_col(node);
      mat &current_nodes_positions = current_pars.positions.slice(t);

      pthread_mutex_lock(lock);
      pthread_cond_wait(position_proposed,lock);
      pthread_mutex_unlock(lock);

      for(n_n node_2 = 0; node_2 < n_nodes; node_2++){
        if(node == node_2) continue;
        if(node_2 % 4 != th_id) continue;
        Q_Q &node_2_type = structure_pars.node_type_each[node_2];

        dist = norm(proposal_position - current_nodes_positions.col(node_2));
        if(node < node_2) proposal_dist.at(node,node_2,t) = dist;
        else proposal_dist.at(node_2,node,t) = dist;

        // node -> node_2
        eta = eta_cal(dist,current_pars.betas.at(node_type,node_2_type),current_pars.radius.at(node_type,node_2_type));
        proposal_llh.at(node,node_2,t) = llh_cal(structure_pars.adj_matrix.at(node,node_2,t),eta);
        // the other way
        eta = eta_cal(dist,current_pars.betas.at(node_2_type,node_type),current_pars.radius.at(node_2_type,node_type));
        proposal_llh.at(node_2,node,t) = llh_cal(structure_pars.adj_matrix.at(node_2,node,t),eta);
      }

      pthread_cond_signal(llh_updated);
    }
  }
  return;
}
*/

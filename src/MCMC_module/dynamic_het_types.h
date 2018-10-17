#ifndef DYNAMIC_HET_TYPES
#define DYNAMIC_HET_TYPES

// n_n meaning non-negative
typedef unsigned int n_n;
// just to make it cute..
typedef char Q_Q;

enum reg_type{
  LOGISTIC, POISSON, TEST, LOGISTIC_F, NOREG
};

enum edge_status{
  // for logistic regression only
  LINKED = 1,
  UNLINKED = 0,

  MISSING = 9,
  MISSING_1 = -91,
  MISSING_0 = -90,

  REMOVED = 10
};

struct network_structure{
  Cube<Q_Q> adj_matrix;
  Mat<Q_Q> node_removed;
  Col<Q_Q> node_type_each;
  Q_Q space_dim;
  Mat<char> regression_type;
  Q_Q derived;
  Q_Q removed_info_derived;

  // derived quantities
  n_n time_length;
  n_n n_nodes;
  Q_Q n_types;
  Col<n_n> n_nodes_each;
  Col<n_n> each_type_start;
  Col<n_n> nodes_rearrange;

  n_n available_nodes_sum;
};

struct network_model_pars{
  double s2;
  double t2;
  cube positions;
  mat betas;
  mat radius;
};

/*
struct model_pars_proposal{
  double s2;
  double t2;
  vec positions;
  mat betas;
  mat radius;
}
*/

struct network_update_pars{
  n_n total_steps;
  n_n current_step;
  n_n burn_in_steps;
  vec positions_step_size;
  double beta_step_size;
  double radius_kappa;
};

struct network_prior_pars{
  double s2_scale;
  double t2_scale;
  double beta_scale;
};

struct acceptance_stat{
  Mat<n_n> positions_accepted;
  Mat<n_n> betas_accepted;
  n_n radius_accepted;
  n_n s2_accepted;
  n_n t2_accepted;
};

struct network_cache{
  cube dist_matrix;
  cube log_likelihood_mat;
};

struct network_info{
  network_structure p_structure;
  network_prior_pars p_prior;
  network_update_pars p_update;
};

struct thread_info{
  void *network_pt;
  Q_Q thread_id;
};

struct thread_info_b{
  void *network_pt;
  Q_Q thread_id;

  Q_Q type_1;
  Q_Q type_2;
};

struct thread_info_p{
  void *network_pt;
  Q_Q thread_id;

  Col<n_n> *update_seq;
  pthread_cond_t *position_proposed;
  pthread_cond_t *llh_updated;
  pthread_mutex_t *lock;
};
#endif

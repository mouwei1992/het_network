#ifndef GEN_POS
#define GEN_POS

class cluster_generator{
private:
  cluster_info generator_info;
  int generated;
  cube positions;
public:
  cluster_generator(cluster_info&);
  void save_positions(const string &) const;
  cube dup_positions(void) const;
  void generate(void);
};

cluster_generator::cluster_generator(cluster_info& cl_info){
  generated = 0;
  generator_info = cl_info;
}

void cluster_generator::generate(void){
  arma_rng::set_seed_random();

  Col<n_n> &n_nodes_each_cluster = generator_info.n_nodes_each_cluster;
  n_n &t_length = generator_info.time_length;
  Q_Q &space_dim = generator_info.space_dimension;

  n_n n_nodes = accu(n_nodes_each_cluster);
  Q_Q n_clusters = n_nodes_each_cluster.n_elem;
  positions = cube(space_dim, n_nodes, t_length);

  n_n current_node = 0;
  for(Q_Q cluster = 0; cluster < n_clusters; cluster++){
    vec cluster_center(space_dim, fill::randn);
    cluster_center *= sqrt(generator_info.var_clusters);
    for(n_n node = 0; node < n_nodes_each_cluster[cluster]; node++){
      positions.slice(0).unsafe_col(current_node).fill(fill::randn);
      positions.slice(0).unsafe_col(current_node) *= sqrt(generator_info.var_points);
      positions.slice(0).unsafe_col(current_node) += cluster_center;

      for(n_n t = 1; t < t_length; t++){
        vec walk(space_dim, fill::randn);
        walk *= generator_info.step_size;
        positions.slice(t).unsafe_col(current_node) = walk + positions.slice(t - 1).unsafe_col(current_node);
      }

      current_node++;
    }
  }

  return;
}

cube cluster_generator::dup_positions(void) const {
  return positions;
}

void cluster_generator::save_positions(const string &file_name_start) const{
  n_n t_length = generator_info.time_length;

  for(n_n t = 0; t < t_length; t++)
    mat_to_csv(file_name_start + to_string(t) + ".csv", positions.slice(t));
}

#endif

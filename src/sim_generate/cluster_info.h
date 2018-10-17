#ifndef CLUSTER_INFO
#define CLUSTER_INFO

struct cluster_info{
  Col<n_n> n_nodes_each_cluster;
  n_n time_length;
  Q_Q space_dimension;

  double var_clusters;
  double var_points;
  double step_size;
};

#endif

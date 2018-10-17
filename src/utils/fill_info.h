#ifndef FILL_INFO
#define FILL_INFO

void fill_structure_info(network_structure &info){
  // info available: adjacency matrix & node type
  if(info.derived) return;

  info.time_length = (info.adj_matrix).n_slices;
  info.n_nodes = (info.adj_matrix).n_cols;

  Col<Q_Q> unique_types = unique(info.node_type_each);
  info.n_types = unique_types.n_elem;

  (info.n_nodes_each) = Col<n_n>(info.n_types);
  (info.each_type_start) = Col<n_n>(info.n_types + 1);
  (info.nodes_rearrange) = Col<n_n>(info.n_nodes);


  // filing removed information
  if(!info.removed_info_derived){
    info.node_removed = Mat<Q_Q> (info.n_nodes, info.time_length);
    info.node_removed.zeros();
    info.removed_info_derived = 1;
  }

  // if(info.regression_type == LOGISTIC){
  //   for(n_n node = 0; node < info.n_nodes; node++){
  //     for(n_n t = 0; t < info.time_length; t++)
  //       info.node_removed.at(node,t) =
  //       (info.adj_matrix.at(node,node,t) == REMOVED)? 1:0;
  //   }
  // }

  // info.available_nodes_sum = info.n_nodes * (info.time_length - 1) - accu(info.node_removed);
  // finding parameter for t^2 update
  info.available_nodes_sum = 0;
  Mat<Q_Q>& rmed = info.node_removed;
  for(n_n node = 0; node < info.n_nodes; node++){
    for(n_n t = 1; t < info.time_length; t++){
      if(!rmed.at(node,t) && !rmed.at(node,t - 1)) info.available_nodes_sum++;
    }
  }


  (info.n_nodes_each).zeros();
  (info.each_type_start).zeros();
  for(n_n node = 0; node < info.n_nodes; node++) ((info.n_nodes_each)((info.node_type_each)(node))) ++;
  for(Q_Q type = 1; type < info.n_types + 1; type++) (info.each_type_start)(type) = (info.each_type_start)(type - 1) + (info.n_nodes_each)(type - 1);

  for(n_n node = 0; node < info.n_nodes; node++){
    Q_Q &node_type = (info.node_type_each)(node);
    n_n &type_start = (info.each_type_start)(node_type);
    (info.nodes_rearrange)(type_start) = node;
    type_start++;
  }

  (info.each_type_start).zeros();
  for(Q_Q type = 1; type < info. n_types + 1; type++) (info.each_type_start)(type) = (info.each_type_start)(type - 1) + (info.n_nodes_each)(type - 1);
  info.derived = 0x01;
  return;
}

#endif

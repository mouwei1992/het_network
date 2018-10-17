#ifndef READ_C_INFO
#define READ_C_INFO

cluster_info read_c_info(const string &info_file_name){
  FILE* fp = fopen(info_file_name.c_str(), "r");
  char readBuffer[65536];
  FileReadStream is(fp, readBuffer, sizeof(readBuffer));

  Document d;
  d.ParseStream(is);

  cluster_info cl_info;

  assert(d.HasMember("n_nodes_each_cluster"));
  const Value& n_nodes_each_cluster = d["n_nodes_each_cluster"];
  assert(n_nodes_each_cluster.IsArray());
  cl_info.n_nodes_each_cluster = Col<n_n>(n_nodes_each_cluster.Size());

  for(Q_Q cl = 0; cl < n_nodes_each_cluster.Size(); cl++) cl_info.n_nodes_each_cluster[cl] = n_nodes_each_cluster[cl].GetInt();

  assert(d.HasMember("time_length"));
  cl_info.time_length = d["time_length"].GetInt();

  assert(d.HasMember("space_dimension"));
  cl_info.space_dimension = d["space_dimension"].GetInt();

  assert(d.HasMember("var_clusters"));
  cl_info.var_clusters = d["var_clusters"].GetDouble();

  assert(d.HasMember("var_points"));
  cl_info.var_points = d["var_points"].GetDouble();

  assert(d.HasMember("step_size"));
  cl_info.step_size = d["step_size"].GetDouble();

  return cl_info;
}

#endif

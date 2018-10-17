#ifndef PARSE_CONFIG
#define PARSE_CONFIG

network_info parse_config(const string &config_file_name, const string& data_path = ""){
  FILE* fp = fopen(config_file_name.c_str(), "r");
  char readBuffer[65536];
  FileReadStream is(fp, readBuffer, sizeof(readBuffer));

  Document d;
  d.ParseStream(is);

  // cannot change const string
  // if(d.HasMember("data_path")) data_path = d["data_path"].GetString();

  network_structure n_str;
  // read adj adj_matrix
  assert(d.HasMember("adj_mat_files"));

  const Value& adj_mat_files = d["adj_mat_files"];

  assert(adj_mat_files.IsArray());
  Mat<Q_Q> adj_mat_static = csv_to_mat<Q_Q>(adj_mat_files[0].GetString(), data_path);
  n_n n_nodes = adj_mat_static.n_cols;
  Cube<Q_Q> adj_mat(n_nodes, n_nodes, adj_mat_files.Size());
  adj_mat.slice(0) = adj_mat_static;

  for(n_n t = 1; t < adj_mat_files.Size(); t++){
    adj_mat_static = csv_to_mat<Q_Q>(adj_mat_files[t].GetString(), data_path);
    adj_mat.slice(t) = adj_mat_static;
  }
  n_str.adj_matrix = adj_mat;

  // read node types
  assert(d.HasMember("types_file"));
  Mat<Q_Q> types = csv_to_mat<Q_Q>(d["types_file"].GetString(), data_path);
  assert(types.is_rowvec() || types.is_colvec());
  if(types.is_rowvec()) types = types.t();
  n_str.node_type_each = types;

  // set other structure
  assert(d.HasMember("space_dim"));
  n_str.space_dim = d["space_dim"].GetInt();

  n_str.derived = 0;

  network_prior_pars n_pp;
  // configure prior parameters
  assert(d.HasMember("s2_scale"));
  n_pp.s2_scale = d["s2_scale"].GetDouble();
  assert(d.HasMember("t2_scale"));
  n_pp.t2_scale = d["t2_scale"].GetDouble();
  assert(d.HasMember("beta_scale"));
  n_pp.beta_scale = d["beta_scale"].GetDouble();

  network_update_pars n_up;
  // configure MCMC parameters
  assert(d.HasMember("total_steps"));
  n_up.total_steps = d["total_steps"].GetInt();
  assert(d.HasMember("burn_in_steps"));
  n_up.burn_in_steps = d["burn_in_steps"].GetInt();

  n_up.current_step = 0;

  assert(d.HasMember("position_step_size"));
  const Value& position_step_size = d["position_step_size"];
  assert(position_step_size.IsArray());
  n_up.positions_step_size = vec(position_step_size.Size());
  for(Q_Q type = 0; type < position_step_size.Size(); type++){
    n_up.positions_step_size[type] = position_step_size[type].GetDouble();
  }

  assert(d.HasMember("beta_step_size"));
  n_up.beta_step_size = d["beta_step_size"].GetDouble();

  assert(d.HasMember("radius_kappa"));
  n_up.radius_kappa = d["radius_kappa"].GetDouble();

  fclose(fp);

  // return the network
  network_info n_info = {n_str, n_pp, n_up};
  return n_info;
}

#endif

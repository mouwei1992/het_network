#ifndef SAVE_RESULT
#define SAVE_RESULT

void print_pars(ostream& ost, const network_model_pars &md_par){
  // does not print out positions
  ost << "model s2: " << md_par.s2 << endl;
  ost << "model t2: " << md_par.t2 << endl;
  ost << "model betas: " << endl << md_par.betas << endl;
  ost << "model radius: " << endl << md_par.radius << endl;

  return;
}

void save_result(const network_model_pars &md_par, const string &par_file_name ="", const string &position_file_name = ""){
  // if not specify file name for saving parameters, simply print on screen
  if(par_file_name.empty()) print_pars(cout, md_par);
  else{
    ofstream file_writer;
    // warning: will overwrite data
    file_writer.open(par_file_name, std::ofstream::out | std::ofstream::app);
    print_pars(file_writer, md_par);
  }

  if(!position_file_name.empty()){
    const cube& positions = md_par.positions;
    for(n_n t = 0; t < positions.n_slices; t++){
      mat_to_csv(position_file_name + to_string(t) + ".csv", positions.slice(t));
    }
  }
  return;
}

#endif

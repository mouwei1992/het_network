#ifndef SAVE_STEPSIZE
#define SAVE_STEPSIZE

void print_stepsize(ostream& ost, const network_update_pars &up_par){
  // does not print out positions
  ost << "total steps: " << up_par.total_steps << endl;
  ost << "current step: " << up_par.current_step << endl;
  ost << "burn-in steps: " << up_par.burn_in_steps << endl;
  ost << "step size for each type of nodes: " << endl << up_par.positions_step_size << endl;
  ost << "beta step size: " << up_par.beta_step_size << endl;
  ost << "radius kappa: " << up_par.radius_kappa << endl;

  return;
}

void save_stepsize(const network_update_pars &up_par, const string &par_file_name =""){
  if(par_file_name.empty()) print_stepsize(cout, up_par);
  else{
    ofstream file_writer;
    // do not overwrite data
    file_writer.open(par_file_name, std::ofstream::out | std::ofstream::app);
    print_stepsize(file_writer, up_par);
  }
}

#endif

#ifndef SAVE_PRIOR
#define SAVE_PRIOR

void print_prior(ostream& ost, const network_prior_pars &pr_par){
  // does not print out positions
  ost << "s2 scale: " << pr_par.s2_scale << endl;
  ost << "t2 scale: " << pr_par.t2_scale << endl;
  ost << "beta scale: " << pr_par.beta_scale << endl;

  return;
}

void save_prior(const network_prior_pars &pr_par, const string &par_file_name =""){
  if(par_file_name.empty()) print_prior(cout, pr_par);
  else{
    ofstream file_writer;
    // do not overwrite data
    file_writer.open(par_file_name, std::ofstream::out | std::ofstream::app);
    print_prior(file_writer, pr_par);
  }
}
#endif

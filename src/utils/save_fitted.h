#ifndef SAVE_FITTED
#define SAVE_FITTED

void save_fitted(const cube &fitted, const string &file_start){
  const n_n t_length = fitted.n_slices;

  for(n_n t = 0; t < t_length; t++){
    mat_to_csv(file_start + to_string(t) + ".csv", fitted.slice(t));
  }

  return;
}

#endif

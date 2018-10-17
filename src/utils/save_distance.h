#ifndef SAVE_DISTANCE
#define SAVE_DISTANCE

void save_distance(const cube &distance, const string &file_start){
  const n_n t_length = distance.n_slices;

  for(n_n t = 0; t < t_length; t++){
    mat_to_csv(file_start + to_string(t) + ".csv", distance.slice(t));
  }

  return;
}

#endif

#ifndef SAVE_ACCEPT
#define SAVE_ACCEPT

void save_acceptance(const dy_het_network &dhn, const string file_start = ""){

  if(file_start.empty()){
    dhn.show_acceptance();
  }else{
    ofstream file_writer;
    // warning: will overwrite data
    file_writer.open(file_start, std::ofstream::out | std::ofstream::trunc);
    dhn.show_acceptance(file_writer);
  }

  return;
}

#endif

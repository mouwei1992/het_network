#ifndef CSV_MAT_H
#define CSV_MAT_H

/*
Assuming the data is stored in square matrix, comma seperated
*/

template<typename T>
Mat<T> csv_to_mat(const string &file_name, const string &path_name = ""){
  ifstream file_reader;
  cout << "opening " << (path_name + file_name) <<endl;
  file_reader.open(path_name + file_name);
  if(!file_reader.good()) cout<<" file " << (path_name + file_name) << " read error! "<<endl;

  string file_line_buf;

  getline(file_reader,file_line_buf);

  regex num_reg("([0-9]+)");
  sregex_iterator current_match(file_line_buf.begin(), file_line_buf.end(), num_reg);
  sregex_iterator last_match;

  n_n n_elem = 0;
  smatch match;

  while(current_match != last_match){
    match = *current_match;
    current_match++;
    n_elem++;
  }

  n_n n_lines = 1;
  while(getline(file_reader,file_line_buf)) n_lines++;
  file_reader.close();

  file_reader.open(path_name + file_name);
  Mat<T> file_in_mat(n_lines,n_elem);

  for(n_n line = 0; line < n_lines; line++){
    getline(file_reader,file_line_buf);
    current_match = sregex_iterator(file_line_buf.begin(), file_line_buf.end(), num_reg);
    for(n_n elem = 0; elem < n_elem; elem++){
      // skip blank line
      if(current_match == last_match){
        cout << "a blank line in csv file" << endl;
        break;
      }

      match = *current_match;
      file_in_mat(line, elem) = stoi(match.str());
      current_match++;
    }
  }

  file_reader.close();

  return file_in_mat;
}

/*
writting data to csv file in the same format
*/
template<typename T>
void mat_to_csv(const string &file_name, Mat<T> file_in_mat){
  ofstream file_writer;
  // warning: will overwrite data
  file_writer.open(file_name, std::ofstream::out | std::ofstream::trunc);

  n_n mat_nrow = file_in_mat.n_rows;
  n_n mat_ncol = file_in_mat.n_cols;

  for(n_n n_row = 0; n_row < mat_nrow; n_row++){
    for(n_n n_col = 0; n_col < mat_ncol; n_col++){
      file_writer << double(file_in_mat(n_row, n_col));
      if(n_col != mat_ncol - 1) file_writer << ',';
    }
    file_writer << '\n';
  }
  file_writer.close();
}

template<typename T>
void mat_to_csv_notrunc(const string &file_name, Mat<T> file_in_mat){
  ofstream file_writer;
  // warning: will overwrite data
  file_writer.open(file_name, std::ofstream::out | std::ofstream::ate);

  n_n mat_nrow = file_in_mat.n_rows;
  n_n mat_ncol = file_in_mat.n_cols;

  for(n_n n_row = 0; n_row < mat_nrow; n_row++){
    for(n_n n_col = 0; n_col < mat_ncol; n_col++){
      file_writer << double(file_in_mat(n_row, n_col));
      if(n_col != mat_ncol - 1) file_writer << ',';
    }
    file_writer << '\n';
  }
  file_writer.close();
}

#endif

#include "../MCMC_module/file_utils.h"

using namespace std;
using namespace arma;

int main(void){
  string CLASSROOM_PATH = "../../data_file/classroom/klas12b-net-";
  mat new_mat(3,3,fill::ones);
  Mat<Q_Q> adj_mat_1 = csv_to_mat<Q_Q>("3.dat", CLASSROOM_PATH);
  mat_to_csv("test_file.csv", adj_mat_1);
  return 0;
}

#include "../utils/dy_het_generator.cpp"
#include "../sim_generate/gen_pos.cpp"

int main(){
  cluster_info c1_info = read_c_info("c1_info.json");
  cluster_generator c1_generator(c1_info);

  c1_generator.generate();

  
}

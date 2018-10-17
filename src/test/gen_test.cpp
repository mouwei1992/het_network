#include "../sim_generate/gen_pos.cpp"

int main(){
  cluster_info c1_info = read_c_info("c1_info.json");
  cluster_generator c1_generator(c1_info);

  c1_generator.generate();

  c1_generator.save_positions("c1_gen");
  return 0;
}

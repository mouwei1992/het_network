#ifndef RANDOM_DISTNS
#define RANDOM_DISTNS


double r_hf_cauchy(double &scale){
  double ru = randu();
  return tan(ru * (datum::pi) / 2.)  * scale;
}

int r_poisson(const double &poisson_mean){
  int rp = 0;
  double ru = randu();
  double sum_exp = -log(ru);

  while(sum_exp < poisson_mean){
    rp++;
    ru = randu();
    sum_exp -= log(ru);
  }

  return rp;
}

template<typename data_storage>
void r_dirichlet(data_storage &concentration){
  // changing input on spot
  concentration.for_each([](double &elem) {elem = elem==0? 0:randg(distr_param(elem,1.));});
  double c_sum = accu(concentration);
  concentration.for_each([&c_sum](double &elem) {elem /= c_sum;});
}

#endif

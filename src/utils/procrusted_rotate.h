#ifndef PROCRUSTES_ROTATE
#define PROCRUSTES_ROTATE\

mat procrusted_rotate(const mat& target_mat, const mat& rot_matrix){
  // assuming the 2 matrices are both centered
  mat u,v;
  vec s;
  svd(u, s, v, target_mat * rot_matrix.t());

  // rotation = u * t(v)
  return u * v.t();
}

#endif

#include "../utils/procrusted_rotate.cpp"

int main(){
  Mat<Q_Q> removed_mat(10, 1, fill::zeros);
  removed_mat(0,0) = 1;
  cout << find(removed_mat.col(0) == 0);
  mat a(2,10, fill::randu);
  a.each_col() -= mean(a, 1);

  mat rot(2,2, fill::zeros);
  cout << rot;
  rot(0,0) = rot(1,1) = cos(0.3);
  rot(0,1) = -sin(0.3);
  rot(1,0) = sin(0.3);

  mat b(a);
  b.each_col([rot](vec &col){col = rot * col;});

  rot = procrusted_rotate(a, b);
  b.each_col([rot](vec &col){col = rot * col;});
  cout << a.cols(find(removed_mat.col(0) == 0));
  return 0;
}

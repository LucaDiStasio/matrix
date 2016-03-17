#include "matrix.h"

//=====================  BODY  ==========================

//-----------------
// Constructors 
//-----------------

template<class T>
matrix<T>::matrix(){}

template<class T>
matrix<T>::matrix(int rows){
  N = rows;
  M = rows;
  A.resize(N);
  for(unsigned int i=0; i<N; i++){
    A[i].resize(M);
    for(unsigned int j=0; j<M; j++){
      A[i][j] = 0;
    }
  }
  measure = false;
}

template<class T>
matrix<T>::matrix(int rows, bool meas){
  N = rows;
  M = rows;
  A.resize(N);
  for(unsigned int i=0; i<N; i++){
    A[i].resize(M);
    for(unsigned int j=0; j<M; j++){
      A[i][j] = 0;
    }
  }
  measure = meas;
}

template<class T>
matrix<T>::matrix(int rows, int columns){
  N = rows;
  M = columns;
  A.resize(N);
  for(unsigned int i=0; i<N; i++){
    A[i].resize(M);
    for(unsigned int j=0; j<M; j++){
      A[i][j] = 0;
    }
  }
  measure = false;
}

template<class T>
matrix<T>::matrix(int rows, int columns, bool meas){
  N = rows;
  M = columns;
  A.resize(N);
  for(unsigned int i=0; i<N; i++){
    A[i].resize(M);
    for(unsigned int j=0; j<M; j++){
      A[i][j] = 0;
    }
  }
  measure = meas;
}

template<class T>
matrix<T>::matrix(T inp[], int rows){
  N = rows;
  M = rows;
  A.resize(N);
  for(unsigned int i=0; i<N; i++){
    A[i].resize(M);
    for(unsigned int j=0; j<M; j++){
      A[i][j] = inp[j + i*N];
    }
  }
  measure = false;
}

template<class T>
matrix<T>::matrix(T inp[], int rows, bool meas){
  N = rows;
  M = rows;
  A.resize(N);
  for(unsigned int i=0; i<N; i++){
    A[i].resize(M);
    for(unsigned int j=0; j<M; j++){
      A[i][j] = inp[j + i*N];
    }
  }
  measure = meas;
}
  
template<class T>
matrix<T>::matrix(T inp[], int rows, int columns){
  N = rows;
  M = columns;
  A.resize(N);
  for(unsigned int i=0; i<N; i++){
    A[i].resize(M);
    for(unsigned int j=0; j<M; j++){
      A[i][j] = inp[j + i*N];
    }
  }
  measure = false;
}
  
template<class T>
matrix<T>::matrix(T inp[], int rows, int columns, bool meas){
  N = rows;
  M = columns;
  A.resize(N);
  for(unsigned int i=0; i<N; i++){
    A[i].resize(M);
    for(unsigned int j=0; j<M; j++){
      A[i][j] = inp[j + i*N];
    }
  }
  measure = meas;
}

template<class T>
matrix<T>::matrix(vector<vector<T> > inp){
  A = inp;
  N = inp.size();
  M = inp[0].size();
  measure = false;
}

template<class T>
matrix<T>::matrix(vector<vector<T> > inp, bool meas){
  A = inp;
  N = inp.size();
  M = inp[0].size();
  measure = meas;
}

//-----------------
//   Destructor 
//-----------------

template<class T>
matrix<T>::~matrix(){}

//-----------------
// General tools 
//-----------------

template<class T>
void matrix<T>::clear(){
  for(unsigned int i=0; i<N; i++){
    A[i].clear();
  }
  A.clear();
}

//-----------------
// Linear Algebra 
//-----------------

template<class T>
void matrix<T>::init_random(double a, double b){
  random_device rd;
  mt19937 eng(rd());
  uniform_real_distribution<> dist(a, b);
  for(unsigned int i=0; i<N; i++){
    for(unsigned int j=0; j<M; j++){
      A[i][j] = dist(eng);
    }
  }
}
  
template<class T>
void matrix<T>::identity(){
  if(!is_square()){
    cout << "The matrix is rectangular. Identity matrix cannot be constructed." << endl;
  }else{
    for(unsigned int i=0; i<N; i++){
      A[i][i] = 1;
    }
  }
}

template<class T>
void matrix<T>::identity(T value){
  if(!is_square()){
    cout << "The matrix is rectangular. Identity matrix cannot be constructed." << endl;
  }else{
    for(unsigned int i=0; i<N; i++){
      A[i][i] = value;
    }
  }
}

template<class T>
void matrix<T>::ones(){
  for(unsigned int i=0; i<N; i++){
    for(unsigned int j=0; j<M; j++){
      A[i][j] = 1;
    }
  }
}

template<class T>
void matrix<T>::ones(T value){
  for(unsigned int i=0; i<N; i++){
    for(unsigned int j=0; j<M; j++){
      A[i][j] = value;
    }
  }
}

template<class T>
void matrix<T>::diag(T vec[], int length, int pos){
  if(!is_square()){
    if(N>M){
      if(pos>0){
        if(length==M-pos){
          for(unsigned int i=0; i<M-pos; i++){
            A[i][i+pos] = vec[i];
          }
        }else{
          cout << "Error: dimension mismatch. Leaving all matrix elements set to zero." << endl;
        }
      }else if(pos<0){
        if(length==M+pos){
          for(unsigned int i=0; i<M+pos; i++){
            A[i-pos][i] = vec[i];
          }
        }else{
          cout << "Error: dimension mismatch. Leaving all matrix elements set to zero." << endl;
        }
      }else{
        if(length==M){
          for(unsigned int i=0; i<M; i++){
            A[i][i] = vec[i];
          }
        }else{
          cout << "Error: dimension mismatch. Leaving all matrix elements set to zero." << endl;
        }
      }
    }else{
      if(pos>0){
        if(length==N-pos){
          for(unsigned int i=0; i<N-pos; i++){
            A[i][i+pos] = vec[i];
          }
        }else{
          cout << "Error: dimension mismatch. Leaving all matrix elements set to zero." << endl;
        }
      }else if(pos<0){
        if(length==N+pos){
          for(unsigned int i=0; i<N+pos; i++){
            A[i-pos][i] = vec[i];
          }
        }else{
          cout << "Error: dimension mismatch. Leaving all matrix elements set to zero." << endl;
        }
      }else{
        if(length==N){
          for(unsigned int i=0; i<N; i++){
            A[i][i] = vec[i];
          }
        }else{
          cout << "Error: dimension mismatch. Leaving all matrix elements set to zero." << endl;
        }
      }
    }
  }else{
    if(pos>0){
      if(length==N-pos){
        for(unsigned int i=0; i<N-pos; i++){
          A[i][i+pos] = vec[i];
        }
      }else{
        cout << "Error: dimension mismatch. Leaving all matrix elements set to zero." << endl;
      }
    }else if(pos<0){
      if(length==N+pos){
        for(unsigned int i=0; i<N+pos; i++){
          A[i-pos][i] = vec[i];
        }
      }else{
        cout << "Error: dimension mismatch. Leaving all matrix elements set to zero." << endl;
      }
    }else{
      if(length==N){
        for(unsigned int i=0; i<N; i++){
          A[i][i] = vec[i];
        }
      }else{
        cout << "Error: dimension mismatch. Leaving all matrix elements set to zero." << endl;
      }
    }
  }
}

template<class T>
void matrix<T>::diag(vector<T> vec, int pos){
  int length = vec.size();
  if(!is_square()){
    if(N>M){
      if(pos>0){
        if(length==M-pos){
          for(unsigned int i=0; i<M-pos; i++){
            A[i][i+pos] = vec[i];
          }
        }else{
          cout << "Error: dimension mismatch. Leaving all matrix elements set to zero." << endl;
        }
      }else if(pos<0){
        if(length==M+pos){
          for(unsigned int i=0; i<M+pos; i++){
            A[i-pos][i] = vec[i];
          }
        }else{
          cout << "Error: dimension mismatch. Leaving all matrix elements set to zero." << endl;
        }
      }else{
        if(length==M){
          for(unsigned int i=0; i<M; i++){
            A[i][i] = vec[i];
          }
        }else{
          cout << "Error: dimension mismatch. Leaving all matrix elements set to zero." << endl;
        }
      }
    }else{
      if(pos>0){
        if(length==N-pos){
          for(unsigned int i=0; i<N-pos; i++){
            A[i][i+pos] = vec[i];
          }
        }else{
          cout << "Error: dimension mismatch. Leaving all matrix elements set to zero." << endl;
        }
      }else if(pos<0){
        if(length==N+pos){
          for(unsigned int i=0; i<N+pos; i++){
            A[i-pos][i] = vec[i];
          }
        }else{
          cout << "Error: dimension mismatch. Leaving all matrix elements set to zero." << endl;
        }
      }else{
        if(length==N){
          for(unsigned int i=0; i<N; i++){
            A[i][i] = vec[i];
          }
        }else{
          cout << "Error: dimension mismatch. Leaving all matrix elements set to zero." << endl;
        }
      }
    }
  }else{
    if(pos>0){
      if(length==N-pos){
        for(unsigned int i=0; i<N-pos; i++){
          A[i][i+pos] = vec[i];
        }
      }else{
        cout << "Error: dimension mismatch. Leaving all matrix elements set to zero." << endl;
      }
    }else if(pos<0){
      if(length==N+pos){
        for(unsigned int i=0; i<N+pos; i++){
          A[i-pos][i] = vec[i];
        }
      }else{
        cout << "Error: dimension mismatch. Leaving all matrix elements set to zero." << endl;
      }
    }else{
      if(length==N){
        for(unsigned int i=0; i<N; i++){
          A[i][i] = vec[i];
        }
      }else{
        cout << "Error: dimension mismatch. Leaving all matrix elements set to zero." << endl;
      }
    }
  }
}

template<class T>
void matrix<T>::hermitian(){
  H.resize(M);
  if(is_same<T,complex<double> >::value){
    for(unsigned int i=0; i<M; i++){
      H[i].resize(N);
      for(unsigned int j=0; j<N; j++){
        H[i][j] = conj(A[j][i]);
      }
    }
  }else{
    for(unsigned int i=0; i<M; i++){
      H[i].resize(N);
      for(unsigned int j=0; j<N; j++){
        H[i][j] = A[j][i];
      }
    }
  }
}

template<class T>
T matrix<T>::trace(){
  T res = 0;
  if(N==M){
    for(unsigned int i=0; i<N; i++){
      res += A[i][i];
    }
  }else{
    cout << "Matrix is rectangular. The method can be applied only to square matrices." << endl;
  }
  return res;
}

template<class T>
T matrix<T>::trace(vector<vector<T> > mat){
  T res = 0;
  if(mat.size()==mat[0].size()){
    for(unsigned int i=0; i<mat.size(); i++){
      res += mat[i][i];
    }
  }else{
    cout << "Matrix is rectangular. The method can be applied only to square matrices." << endl;
  }
  return res;
}

template<class T>
T matrix<T>::trace(matrix<T> mat){
  T res = 0;
  vector<int> sizes = mat.size();
  if(sizes[0]==sizes[1]){
    for(unsigned int i=0; i<sizes[0]; i++){
      res += mat(i,i);
    }
  }else{
    cout << "Matrix is rectangular. The method can be applied only to square matrices." << endl;
  }
  return res;
}

template<class T>
T matrix<T>::matnorm1(){
  T res = 0;
  vector<T> sums;
  sums.resize(M);
  for(unsigned int j=0; j<M; j++){
    sums[j] = 0;
    for(unsigned int i=0; i<N; i++){
      sums[j] += abs(A[i][j]);
    }
  }
  res = sums[0];
  for(unsigned int j=0; j<M; j++){
    if(sums[j]>res){
      res = sums[j];
    }
  }
  return res;
}

template<class T>
T matrix<T>::matnormInf(){
  T res = 0;
  vector<T> sums;
  sums.resize(N);
  for(unsigned int i=0; i<N; i++){
    sums[i] = 0;
    for(unsigned int j=0; j<M; j++){
      sums[i] += abs(A[i][j]);
    }
  }
  res = sums[0];
  for(unsigned int i=0; i<N; i++){
    if(sums[i]>res){
      res = sums[i];
    }
  }
  return res;
}

template<class T>
T matrix<T>::matnormFrob(){
  hermitian();
  T res = sqrt(trace(prod(A,H)));
  return res;
}

template<class T>
T matrix<T>::vecnorm2(){
  T res = 0;
  for(unsigned int i=0; i<N; i++){
    for(unsigned int j=0; j<M; j++){
      res += pow(A[i][j],2);
    }
  }
  return sqrt(res);
}

template<class T>
vector<vector<T> > matrix<T>::absv(){
  vector<vector<T> > res;
  res.resize(N);
  for(unsigned int i=0; i<N; i++){
    res[i].resize(N);
    for(unsigned int j=0; j<M; j++){
      res[i][j] = abs(A[i][j]);
    }
  }
  return res;
}

template<class T>
vector<vector<T> > matrix<T>::absv(vector<vector<T> > mat){
  int n = mat.size();
  int m = mat[0].size();
  vector<vector<T> > res;
  res.resize(n);
  for(unsigned int i=0; i<n; i++){
    res[i].resize(m);
    for(unsigned int j=0; j<m; j++){
      res[i][j] = abs(mat[i][j]);
    }
  }
  return res;
}

template<class T>
vector<vector<T> > matrix<T>::absv(matrix<T> mat){
  vector<int> sizes = mat.sizes();
  int n = sizes[0];
  int m = sizes[1];
  vector<vector<T> > res;
  res.resize(n);
  for(unsigned int i=0; i<n; i++){
    res[i].resize(m);
    for(unsigned int j=0; j<m; j++){
      res[i][j] = abs(mat(i,j));
    }
  }
  return res;
}

template<class T>
matrix<T> matrix<T>::absm(){
  matrix<T> res(N,M);
  for(unsigned int i=0; i<N; i++){
    for(unsigned int j=0; j<M; j++){
      res(i,j) = abs(A[i][j]);
    }
  }
  return res;
}

template<class T>
matrix<T> matrix<T>::absm(vector<vector<T> > mat){
  int n = mat.size();
  int m = mat[0].size();
  matrix<T> res(n,m);
  for(unsigned int i=0; i<n; i++){
    for(unsigned int j=0; j<m; j++){
      res(i,j) = abs(mat[i][j]);
    }
  }
  return res;
}

template<class T>
matrix<T> matrix<T>::absm(matrix<T> mat){
  vector<int> sizes = mat.sizes();
  int n = sizes[0];
  int m = sizes[1];
  matrix<T> res(n,m);
  for(unsigned int i=0; i<n; i++){
    for(unsigned int j=0; j<m; j++){
      res(i,j) = abs(mat(i,j));
    }
  }
  return res;
}

template<class T>
vector<int> matrix<T>::min(){
  vector<int> res;
  res.resize(2);
  T value = A[0][0];
  res[0] = 0;
  res[1] = 0;
  for(unsigned int i=0; i<N; i++){
    for(unsigned int j=0; j<M; j++){
      if(A[i][j]<value){
        value = A[i][j];
        res[0] = i;
        res[1] = j;
      }
    }
  }
  return res;
}

template<class T>
vector<int> matrix<T>::min(vector<vector<T> > mat){
  int R = mat.size();
  int Q = mat[0].size();
  vector<int> res;
  res.resize(2);
  T value = mat[0][0];
  res[0] = 0;
  res[1] = 0;
  for(unsigned int i=0; i<R; i++){
    for(unsigned int j=0; j<Q; j++){
      if(mat[i][j]<value){
        value = mat[i][j];
        res[0] = i;
        res[1] = j;
      }
    }
  }
  return res;
}
    
template<class T>
vector<int> matrix<T>::min(matrix<T> mat){
  vector<int> sizes = mat.sizes();
  vector<int> res;
  res.resize(2);
  T value = mat(0,0);
  res[0] = 0;
  res[1] = 0;
  for(unsigned int i=0; i<sizes[0]; i++){
    for(unsigned int j=0; j<sizes[1]; j++){
      if(mat(i,j)<value){
        value = mat(i,j);
        res[0] = i;
        res[1] = j;
      }
    }
  }
  return res;
}

template<class T>
vector<int> matrix<T>::max(){
  vector<int> res;
  res.resize(2);
  T value = A[0][0];
  res[0] = 0;
  res[1] = 0;
  for(unsigned int i=0; i<N; i++){
    for(unsigned int j=0; j<M; j++){
      if(A[i][j]>value){
        value = A[i][j];
        res[0] = i;
        res[1] = j;
      }
    }
  }
  return res;
}

template<class T>
vector<int> matrix<T>::max(vector<vector<T> > mat){
  int R = mat.size();
  int Q = mat[0].size();
  vector<int> res;
  res.resize(2);
  T value = mat[0][0];
  res[0] = 0;
  res[1] = 0;
  for(unsigned int i=0; i<R; i++){
    for(unsigned int j=0; j<Q; j++){
      if(mat[i][j]>value){
        value = mat[i][j];
        res[0] = i;
        res[1] = j;
      }
    }
  }
  return res;
}

template<class T>
vector<int> matrix<T>::max(matrix<T> mat){
  vector<int> sizes = mat.sizes();
  vector<int> res;
  res.resize(2);
  T value = mat(0,0);
  res[0] = 0;
  res[1] = 0;
  for(unsigned int i=0; i<sizes[0]; i++){
    for(unsigned int j=0; j<sizes[1]; j++){
      if(mat(i,j)>value){
        value = mat(i,j);
        res[0] = i;
        res[1] = j;
      }
    }
  }
  return res;
}

template<class T>
matrix<T> matrix<T>::vec_min(int flag = 0){
  if(flag){
    matrix<T> res(N,1);
    for(unsigned int i=0; i<N; i++){
      res(i) = A[i][0];
      for(unsigned int j=0; j<M; j++){
        if(A[i][j]<res(i)){
          res(i) = A[i][j];
        }
      }
    }
    return res;
  }else{
    matrix<T> res(1,M);
    for(unsigned int j=0; j<M; j++){
      res(j) = A[0][j];
      for(unsigned int i=0; i<N; i++){
        if(A[i][j]<res(j)){
          res(j) = A[i][j];
        }
      }
    }
    return res;
  }
}

template<class T>
matrix<T> matrix<T>::vec_min(vector<vector<T> > mat, int flag = 0){
  int R = mat.size();
  int Q = mat[0].size();
  if(flag){
    matrix<T> res(R,1);
    for(unsigned int i=0; i<R; i++){
      res(i) = mat[i][0];
      for(unsigned int j=0; j<Q; j++){
        if(mat[i][j]<res(i)){
          res(i) = mat[i][j];
        }
      }
    }
    return res;
  }else{
    matrix<T> res(1,Q);
    for(unsigned int j=0; j<R; j++){
      res(j) = mat[0][j];
      for(unsigned int i=0; i<Q; i++){
        if(mat[i][j]<res(j)){
          res(j) = mat[i][j];
        }
      }
    }
    return res;
  }
}

template<class T>
matrix<T> matrix<T>::vec_min(matrix<T> mat, int flag = 0){
  vector<int> sizes = mat.sizes();
  if(flag){
    matrix<T> res(sizes[0],1);
    for(unsigned int i=0; i<sizes[0]; i++){
      res(i) = mat(i,0);
      for(unsigned int j=0; j<sizes[1]; j++){
        if(mat(i,j)<res(i)){
          res(i) = mat(i,j);
        }
      }
    }
    return res;
  }else{
    matrix<T> res(1,sizes[1]);
    for(unsigned int j=0; j<sizes[1]; j++){
      res(j) = mat(0,j);
      for(unsigned int i=0; i<sizes[0]; i++){
        if(mat(i,j)<res(j)){
          res(j) = mat(i,j);
        }
      }
    }
    return res;
  }
}

template<class T>
matrix<T> matrix<T>::vec_max(int flag = 0){
  if(flag){
    matrix<T> res(N,1);
    for(unsigned int i=0; i<N; i++){
      res(i) = A[i][0];
      for(unsigned int j=0; j<M; j++){
        if(A[i][j]>res(i)){
          res(i) = A[i][j];
        }
      }
    }
    return res;
  }else{
    matrix<T> res(1,M);
    for(unsigned int j=0; j<M; j++){
      res(j) = A[0][j];
      for(unsigned int i=0; i<N; i++){
        if(A[i][j]>res(j)){
          res(j) = A[i][j];
        }
      }
    }
    return res;
  }
}

template<class T>
matrix<T> matrix<T>::vec_max(vector<vector<T> > mat, int flag = 0){
  int R = mat.size();
  int Q = mat[0].size();
  if(flag){
    matrix<T> res(R,1);
    for(unsigned int i=0; i<R; i++){
      res(i) = mat[i][0];
      for(unsigned int j=0; j<Q; j++){
        if(mat[i][j]>res(i)){
          res(i) = mat[i][j];
        }
      }
    }
    return res;
  }else{
    matrix<T> res(1,Q);
    for(unsigned int j=0; j<R; j++){
      res(j) = mat[0][j];
      for(unsigned int i=0; i<Q; i++){
        if(mat[i][j]>res(j)){
          res(j) = mat[i][j];
        }
      }
    }
    return res;
  }
}

template<class T>
matrix<T> matrix<T>::vec_max(matrix<T> mat, int flag = 0){
  vector<int> sizes = mat.sizes();
  if(flag){
    matrix<T> res(sizes[0],1);
    for(unsigned int i=0; i<sizes[0]; i++){
      res(i) = mat(i,0);
      for(unsigned int j=0; j<sizes[1]; j++){
        if(mat(i,j)>res(i)){
          res(i) = mat(i,j);
        }
      }
    }
    return res;
  }else{
    matrix<T> res(1,sizes[1]);
    for(unsigned int j=0; j<sizes[1]; j++){
      res(j) = mat(0,j);
      for(unsigned int i=0; i<sizes[0]; i++){
        if(mat(i,j)>res(j)){
          res(j) = mat(i,j);
        }
      }
    }
    return res;
  }
}

template<class T>
matrix<T> matrix<T>::right_prod(matrix<T> mat){
  vector<int> dim = mat.size();
  int P = dim[0];
  int Q = dim[1];
  matrix<T> result(N,Q);
  if(M==P){
    for(unsigned int n=0; n<N; n++){
      for(unsigned int q=0; q<Q; q++){
        result(n,q) = 0;
        for(unsigned int m=0; m<M; m++){
          result(n,q) += A[n][m]*mat(m,q);
        }
      }
    }
  }else{
    cout << "Error: dimension mismatch. Result left to zero.";
  }
  return result;
}

template<class T>
matrix<T> matrix<T>::left_prod(matrix<T> mat){
  vector<int> dim = mat.size();
  int P = dim[0];
  int Q = dim[1];
  matrix<T> result(P,M);
  if(Q==N){
    for(unsigned int p=0; p<P; p++){
      for(unsigned int m=0; m<M; m++){
        result(p,m) = 0;
        for(unsigned int n=0; n<N; n++){
          result(p,m) += mat(p,n)*A[n][m];
        }
      }
    }
  }else{
    cout << "Error: dimension mismatch. Result left to zero.";
  }
  return result;
}

template<class T>
matrix<T> matrix<T>::prod(matrix<T> mat1, matrix<T> mat2){
  vector<int> dim = mat1.size();
  int P = dim[0];
  int Q = dim[1];
  dim.clear();
  dim = mat2.size();
  int R = dim[0];
  int S = dim[1];
  matrix<T> result(P,S);
  if(Q==R){
    for(unsigned int p=0; p<P; p++){
      for(unsigned int s=0; s<S; s++){
        result(p,s) = 0;
        for(unsigned int q=0; q<Q; q++){
          result(p,s) += mat1(p,q)*mat2(q,s);
        }
      }
    }
  }else{
    cout << "Error: dimension mismatch. Result left to zero.";
  }
  return result;
}

template<class T>
vector<vector<T> > matrix<T>::prod(vector<vector<T> > mat1, vector<vector<T> > mat2){
  int P = mat1.size();
  int Q = mat1[0].size();
  int R = mat2.size();
  int S = mat2[0].size();
  vector<vector<T> > result;
  if(Q==R){
    result.resize(P);
    for(unsigned int p=0; p<P; p++){
      result[p].resize(S);
      for(unsigned int s=0; s<S; s++){
        result[p][s] = 0;
        for(unsigned int q=0; q<Q; q++){
          result[p][s] += mat1[p][q]*mat2[q][s];
        }
      }
    }
  }else{
    cout << "Error: dimension mismatch. Result left to zero.";
  }
  return result;
}

template<class T>
T matrix<T>::right_scalar_prod(matrix<T> mat){
  vector<int> dim = mat.size();
  int P = dim[0];
  int Q = dim[1];
  T result = 0;
  if(N==1 && M==P && Q==1){
    for(unsigned int m=0; m<M; m++){
      result += A[0][m]*mat(m,0);
    }
  }else{
    cout << "Error: dimension mismatch. Result left to zero.";
  }
  return result;
}

template<class T>
T matrix<T>::left_scalar_prod(matrix<T> mat){
  vector<int> dim = mat.size();
  int P = dim[0];
  int Q = dim[1];
  T result = 0;
  if(M==1 && N==Q && P==1){
    for(unsigned int n=0; n<N; n++){
      result += mat(0,n)*A[n][0];
    }
  }else{
    cout << "Error: dimension mismatch. Result left to zero.";
  }
  return result;
}

template<class T>
T matrix<T>::scalar_prod(matrix<T> mat1, matrix<T> mat2){
  vector<int> dim = mat1.size();
  int P = dim[0];
  int Q = dim[1];
  dim.clear();
  dim = mat2.size();
  int R = dim[0];
  int S = dim[1];
  T result = 0;
  if(P==1 && Q==R && S==1){
    for(unsigned int q=0; q<Q; q++){
      result += mat1(0,q)*mat1(q,0);
    }
  }else{
    cout << "Error: dimension mismatch. Result left to zero.";
  }
  return result;
}

template<class T>
matrix<T> matrix<T>::right_tensor_prod(matrix<T> mat){
  vector<int> dim = mat.size();
  int R = dim[0];
  int S = dim[1];
  matrix<T> result(N*R,M*S);
  for(unsigned int n=0; n<N; n++){
    for(unsigned int r=0; r<R; r++){
      for(unsigned int m=0; m<M; m++){
        for(unsigned int s=0; s<S; s++){
          result(r+n*R,s+m*S) = A[n][m]*mat(r,s);  
        }
      }  
    }
  }
  return result;
}

template<class T>
matrix<T> matrix<T>::left_tensor_prod(matrix<T> mat){
  vector<int> dim = mat.size();
  int P = dim[0];
  int Q = dim[1];
  matrix<T> result(P*N,Q*M);
  for(unsigned int p=0; p<P; p++){
    for(unsigned int n=0; n<N; n++){
      for(unsigned int q=0; q<Q; q++){
        for(unsigned int m=0; m<M; m++){
          result(n+p*N,m+q*M) = mat(p,q)*A[n][m];  
        }
      }  
    }
  }
  return result;
}

template<class T>
matrix<T> matrix<T>::tensor_prod(matrix<T> mat1, matrix<T> mat2){
  vector<int> dim = mat1.size();
  int P = dim[0];
  int Q = dim[1];
  dim.clear();
  dim = mat2.size();
  int R = dim[0];
  int S = dim[1];
  matrix<T> result(P*R,Q*S);
  for(unsigned int p=0; p<P; p++){
    for(unsigned int r=0; r<R; r++){
      for(unsigned int q=0; q<Q; q++){
        for(unsigned int s=0; s<S; s++){
          result(r+p*R,s+q*S) = mat1(p,q)*mat2(r,s);  
        }
      }  
    }
  }
  return result;
}

template<class T>
vector<vector<T> > matrix<T>::get_cofactor(vector<vector<T> > B, int n, int row, int col){
  int colCount;
  int rowCount;
  vector<vector<T> > C;
  C.resize(n-1);
  for(unsigned int i=0; i<n-1; i++){
    C[i].resize(n-1);
  }
  rowCount = 0;
  for(unsigned int i=0; i<n; i++){
    colCount = 0;
    if(i!=row){
      for(unsigned int j=0; j<n; j++){
        if(j!=col){
          C[rowCount][colCount] = B[i][j];
          colCount++;
        }
      }
      rowCount++;
    }
  }
  return C;
}

template<class T>
matrix<T> matrix<T>::forw_subs(matrix<T> b){
  matrix<T> x(N,1);
  if(is_square(L) && !is_zero(L)){
    if(abs(L[0][0])!=0){
      x(0) = b(0)/L[0][0];
      T sum;
      for(unsigned int i=1; i<N; i++){
        if(abs(L[i][i])!=0){
          sum = 0;
          for(unsigned int j=0; j<i; j++){
            sum += L[i][j]*x(j);
          }
          x(i) = (b(i)-sum)/L[i][i];
        }else{
          cout << "Error: null diagonal element (" << i << ", " << i << ") of L. Process interrupted." << endl;
          return x;
        }
      }
    }else{
        cout << "Error: null diagonal element (" << 0 << ", " << 0 << ") of L. Process interrupted." << endl;
    }
  }else{
    cout << "Error: lower triangular matrix L not square or equal to null matrix. Result left to zero." << endl;
  }
  return x;
}

template<class T>
matrix<T> matrix<T>::forw_subs(matrix<T> Linp, matrix<T> b){
  vector<int> sizes = Linp.size();
  matrix<T> x(sizes[0],1);
  if(is_square(Linp) && !is_zero(Linp)){
    if(abs(Linp(0,0))!=0){
      x(0) = b(0)/Linp(0,0);
      T sum;
      for(unsigned int i=1; i<sizes[0]; i++){
        if(abs(Linp(i,i))!=0){
          sum = 0;
          for(unsigned int j=0; j<i; j++){
            sum += Linp(i,j)*x(j);
          }
          x(i) = (b(i)-sum)/Linp(i,i);
        }else{
          cout << "Error: null diagonal element (" << i << ", " << i << ") of L. Process interrupted." << endl;
          return x;
        }
      }
    }else{
        cout << "Error: null diagonal element (" << 0 << ", " << 0 << ") of L. Process interrupted." << endl;
    }
  }else{
    cout << "Error: lower triangular matrix L not square or equal to null matrix. Result left to zero." << endl;
  }
  return x;
}

template<class T>
vector<T> matrix<T>::forw_subs(vector<T> b){
  vector<T> x;
  x.resize(N);
  if(is_square(L) && !is_zero(L)){
    if(abs(L[0][0])!=0){
      x[0] = b[0]/L[0][0];
      T sum;
      for(unsigned int i=1; i<N; i++){
        if(abs(L[i][i])!=0){
          sum = 0;
          for(unsigned int j=0; j<i; j++){
            sum += L[i][j]*x[j];
          }
          x[i] = (b[i]-sum)/L[i][i];
        }else{
          cout << "Error: null diagonal element (" << i << ", " << i << ") of L. Process interrupted." << endl;
          return x;
        }
      }
    }else{
        cout << "Error: null diagonal element (" << 0 << ", " << 0 << ") of L. Process interrupted." << endl;
    }
  }else{
    cout << "Error: lower triangular matrix L not square or equal to null matrix. Result left to zero." << endl;
  }
  return x;
}

template<class T>
vector<T> matrix<T>::forw_subs(vector<vector<T> > Linp, vector<T> b){
  int n = Linp.size();
  vector<T> x;
  x.resize(n);
  if(is_square(Linp) && !is_zero(Linp)){
    if(abs(Linp[0][0])!=0){
      x[0] = b[0]/Linp[0][0];
      T sum;
      for(unsigned int i=1; i<n; i++){
        if(abs(Linp[i][i])!=0){
          sum = 0;
          for(unsigned int j=0; j<i; j++){
            sum += Linp[i][j]*x[j];
          }
          x[i] = (b[i]-sum)/Linp[i][i];
        }else{
          cout << "Error: null diagonal element (" << i << ", " << i << ") of L. Process interrupted." << endl;
          return x;
        }
      }
    }else{
        cout << "Error: null diagonal element (" << 0 << ", " << 0 << ") of L. Process interrupted." << endl;
    }
  }else{
    cout << "Error: lower triangular matrix L not square or equal to null matrix. Result left to zero." << endl;
  }
  return x;
}

template<class T>
matrix<T> matrix<T>::back_subs(matrix<T> b){
  matrix<T> x(N,1);
  if(is_square(U) && !is_zero(U)){
    if(abs(U[N-1][N-1])!=0){
      x(N-1) = b(N-1)/U[N-1][N-1];
      T sum;
      for(unsigned int i=N-2; i>=0; i--){
        if(abs(U[i][i])!=0){
          sum = 0;
          for(unsigned int j=i+1; j<N; j++){
            sum += U[i][j]*x(j);
          }
          x(i) = (b(i)-sum)/U[i][i];
        }else{
          cout << "Error: null diagonal element (" << i << ", " << i << ") of U. Process interrupted." << endl;
          return x;
        }
      }
    }else{
        cout << "Error: null diagonal element (" << 0 << ", " << 0 << ") of U. Process interrupted." << endl;
    }
  }else{
    cout << "Error: lower triangular matrix U not square or equal to null matrix. Result left to zero." << endl;
  }
  return x;
}

template<class T>
matrix<T> matrix<T>::back_subs(matrix<T> Uinp, matrix<T> b){
  vector<int> sizes = Uinp.size();
  matrix<T> x(sizes[0],1);
  if(is_square(Uinp) && !is_zero(U)){
    if(abs(Uinp(sizes[0]-1,sizes[0]-1))!=0){
      x(sizes[0]-1) = b(sizes[0]-1)/Uinp(sizes[0]-1,sizes[0]-1);
      T sum;
      for(unsigned int i=sizes[0]-2; i>=0; i--){
        if(abs(Uinp(i,i))!=0){
          sum = 0;
          for(unsigned int j=i+1; j<sizes[0]; j++){
            sum += Uinp(i,j)*x(j);
          }
          x(i) = (b(i)-sum)/Uinp(i,i);
        }else{
          cout << "Error: null diagonal element (" << i << ", " << i << ") of U. Process interrupted." << endl;
          return x;
        }
      }
    }else{
        cout << "Error: null diagonal element (" << 0 << ", " << 0 << ") of U. Process interrupted." << endl;
    }
  }else{
    cout << "Error: lower triangular matrix U not square or equal to null matrix. Result left to zero." << endl;
  }
  return x;
}

template<class T>
vector<T> matrix<T>::back_subs(vector<T> b){
  vector<T> x;
  x.resize(N);
  if(is_square(U) && !is_zero(U)){
    if(abs(U[N-1][N-1])!=0){
      x[N-1] = b[N-1]/U[N-1][N-1];
      T sum;
      for(unsigned int i=N-2; i>=0; i--){
        if(abs(U[i][i])!=0){
          sum = 0;
          for(unsigned int j=i+1; j<N; j++){
            sum += U[i][j]*x[j];
          }
          x[i] = (b[i]-sum)/U[i][i];
        }else{
          cout << "Error: null diagonal element (" << i << ", " << i << ") of U. Process interrupted." << endl;
          return x;
        }
      }
    }else{
        cout << "Error: null diagonal element (" << 0 << ", " << 0 << ") of U. Process interrupted." << endl;
    }
  }else{
    cout << "Error: lower triangular matrix U not square or equal to null matrix. Result left to zero." << endl;
  }
  return x;
}

template<class T>
vector<T> matrix<T>::back_subs(vector<vector<T> > Uinp, vector<T> b){
  int n = Uinp.size();
  vector<T> x;
  x.resize(n);
  if(is_square(Uinp) && !is_zero(U)){
    if(abs(Uinp[n-1][n-1])!=0){
      x[n-1] = b[n-1]/Uinp[n-1][n-1];
      T sum;
      for(unsigned int i=n-2; i>=0; i--){
        if(abs(Uinp[i][i])!=0){
          sum = 0;
          for(unsigned int j=i+1; j<n; j++){
            sum += Uinp[i][j]*x[j];
          }
          x[i] = (b[i]-sum)/Uinp[i][i];
        }else{
          cout << "Error: null diagonal element (" << i << ", " << i << ") of U. Process interrupted." << endl;
          return x;
        }
      }
    }else{
        cout << "Error: null diagonal element (" << 0 << ", " << 0 << ") of U. Process interrupted." << endl;
    }
  }else{
    cout << "Error: lower triangular matrix U not square or equal to null matrix. Result left to zero." << endl;
  }
  return x;
}

template<class T>
void matrix<T>::LU(){
  if(is_square()){
    L.resize(N);
    for(unsigned int i=0; i<N; i++){
      L[i].resize(N);
    }
    U.resize(N);
    for(unsigned int i=0; i<N; i++){
      U[i].resize(N);
    }
    T sum;
    for(unsigned int k=0; k<N; k++){
      if(A[k][k]!=0){
        L[k][k] = 1;
        for(unsigned int j=k; j<N; j++){
          sum = 0;
          for(unsigned int r=0; r<k; r++){
            sum += L[k][r]*U[r][j];
          }
          U[k][j] = A[k][j] - sum;
        }
        for(unsigned int i=k+1; i<N; i++){
          sum = 0;
          for(unsigned int r=0; r<k; r++){
            sum += L[i][r]*U[r][k];
          }
          L[i][k] = (A[i][k] - sum)/U[k][k];
        }
      }else{
        cout << "Null diagonal element (" << k << ", " << k << "). Process interrupted." << endl;
        break;
      }
    }
  }else{
    cout << "Matrix is rectangular. The method can be applied only to square matrices." << endl;
  }
}

template<class T>
void matrix<T>::LU_pp(){
  if(is_square()){
    L.resize(N);
    for(unsigned int i=0; i<N; i++){
      L[i].resize(N);
    }
    U.resize(N);
    for(unsigned int i=0; i<N; i++){
      U[i].resize(N);
    }
    T sum;
    U = A;
    identity();
    L = A;
    P = A;
    A = U;
    vector<vector<T> > submatrix;
    vector<int> pivot;
    T rowEl;
    int rowPivot;
    for(unsigned int k=0; k<N-1; k++){
      submatrix.clear();
      submatrix.resize(N-k);
      for(unsigned int i=k; i<N; i++){
        submatrix[i-k].resize(N-k);
        for(unsigned int j=k; j<N; j++){
          submatrix[i-k][j-k] = U[i][j];
        }
      }
      pivot = max(submatrix);
      rowPivot = pivot[0] + k;
      if(rowPivot!=k){
        for(unsigned int i=k; i<N; k++){
          rowEl = 0;
          rowEl = U[k][i];
          U[k][i] = U[rowPivot][i];
          U[rowPivot][i] = rowEl;
        }
        for(unsigned int i=0; i<k; k++){
          rowEl = 0;
          rowEl = L[k][i];
          L[k][i] = L[rowPivot][i];
          L[rowPivot][i] = rowEl;
        }
        for(unsigned int i=0; i<N; k++){
          rowEl = 0;
          rowEl = P[k][i];
          P[k][i] = P[rowPivot][i];
          P[rowPivot][i] = rowEl;
        }
      }
      for(unsigned int j=k+1; j<N; j++){
        L[j][k] = U[j][k]/U[k][k];
        for(unsigned int i=k; i<N; i++){
          U[j][i] += L[j][k]*U[k][i];
        }
      }
    }
  }else{
    cout << "Matrix is rectangular. The method can be applied only to square matrices." << endl;
  }
}

template<class T>
void matrix<T>::LU_fp(){
  if(is_square()){
    T max;
    int ipiv, jpiv;
    vector<vector<T> > B = A;
    A.clear();
    identity();
    vector<vector<T> > I = A;
    P = I;
    A.clear();
    identity(2);
    vector<vector<T> > I2 = A;
    A = B;
    Q = I;
    vector<vector<T> > Minv = I;
    vector<vector<T> > Pk;
    vector<vector<T> > Qk;
    vector<vector<T> > Mk;
    vector<vector<T> > Mkinv;
    vector<vector<T> > submatrix;
    vector<int> pij;
    for(unsigned int k=0; k<N-1; k++){
      submatrix.clear();
      submatrix.resize(N-k);
      for(unsigned int i=k; i<N; i++){
        submatrix[i].resize(N-k);
        for(unsigned int j=k; j<M; j++){
          submatrix[i-k][j-k] = B[i][j];
        }
      }
      pij.clear();
      pij = max(absv(submatrix));
      pij[0] += k;
      pij[1] += k;
      Pk.clear();
      Pk = I;
      Pk(pij[0],pij[0]) = 0;
      Pk(k,k) = 0;
      Pk(k,pij[0]) = 1;
      Pk(pij[0],k) = 1;
      Qk.clear();
      Qk = I;
      Qk(pij[1],pij[1]) = 0;
      Qk(k,k) = 0;
      Qk(k,pij[1]) = 1;
      Qk(pij[1],k) = 1;
      B = prod(Pk,prod(B,Qk));
      Mk.clear();
      Mk = I;
      for(unsigned int i=k+1; i<N; i++){
        Mk[i][k] = -A[i][k]/A[k][k];
      }
      Mkinv = I2;
      for(unsigned int i=0; i<N; i++){
        for(unsigned int j=0; j<M; j++){
          Mkinv[i][j] -= Mk[i][j]; 
        }
      }
      B = prod(Mk,B);
      P = prod(Pk,P);
      Q = prod(Q,Qk);
      Minv = prod(Minv,prod(Pk,Mkinv));
    }
    U.resize(N);
    for(unsigned int i=0; i<N; i++){
      U[i].resize(M);
      for(unsigned int j=i; j<M; j++){
        U[i][j] = B[i][j];
      }
    }
    L = prod(P,Minv);
  }else{
    cout << "Matrix is rectangular. The method can be applied only to square matrices." << endl;
  }
}

template<class T>
void matrix<T>::Cholesky(){
  if(is_square()){
    if(is_symm()){
      if(is_posdef()){
        U.resize(N);
        for(unsigned int i=0; i<N; i++){
          U[i].resize[M];
        }
        T sum;
        U[0][0] = sqrt(A[0][0]);
        for(unsigned int i=1; i<N; i++){
          if(abs(A[i][i])<=0){
            cout << "Null or negative diagonal element (" << i << ", " << i << "). Factorization cannot be applied. Process interrupted." << endl;
            exit(0);
          }else{
            for(unsigned int j=0; j<i; j++){
              sum = 0;
              for(unsigned int k=0; k<j; k++){
                sum += U[i][k]*U[j][k];
              }
              U[i][j] = (A[i][j]-sum)/U[j][j];
            }
            sum = 0;
            for(unsigned int k=0; k<i; k++){
              sum += U[i][k]*U[i][k];
            }
            U[i][i] = sqrt(A[i][i]-sum);
          }
        }
        L.resize(N);
        for(unsigned int i=0; i<N; i++){
          L[i].resize[M];
          for(unsigned int j=0; j<=i; j++){
            L[i][j] = U[j][i];
          }
        }
      }else{
        cout << "Matrix not positive definite. Cholesky factorization cannot be applied." << endl;
      }
    }else{
      cout << "Matrix not symmetric. Cholesky factorization cannot be applied." << endl;
    }
    
  }else{
    cout << "Matrix is rectangular. The method can be applied only to square matrices." << endl;
  }
}

template<class T>
T matrix<T>::recursive_det(vector<vector<T> > B, int n){
  T D = 0;
  unsigned int l;
  if (n<1){
    cout << "Matrix dimension is negative." << endl;
    exit(0);
  }else if (n == 1) {
    D = B[0][0];
  }else if (n == 2) {
    D = B[0][0] * B[1][1] - B[1][0] * B[0][1];
  }else {
    D = 0;
    int i = 0;
    for (int j=0; j<n; j++) {
      D += pow(-1.0,(i+j)%2) * B[i][j] * recursive_det(get_cofactor(B,n,i,j),n-1);
    }
  }
  return D; 
}

template<class T>
T matrix<T>::LU_det(){
  LU_pp();
  int t = 0;
  T detU, detL, detP;
  for(unsigned int i=0; i<N; i++){
    if(abs(P[i][i])==0){
      t++;
    }
  }
  if(t>0){
    detP = pow(-1,t-1);
  }else{
    detP = 1;
  }
  detU = 1;
  detL = 1;
  for(unsigned int i=0; i<N; i++){
    detU *= U[i][i];
    detL *= L[i][i];
  }
  T deter = detP*detU*detL;
  return deter;
}

template<class T>  
vector<vector<T> > matrix<T>::recursive_inv(vector<vector<T> > B, int dim){
  vector<vector<T> > inv_mat;
  if (dim<1){
    cout << "Matrix dimension is negative." << endl;
    exit(0);
  }else if (dim == 1) {
    inv_mat.resize(dim);
    inv_mat[0].resize(dim);
    inv_mat[0][0] = B[0][0];
  }else if (dim == 2) {
    inv_mat.resize(dim);
    inv_mat[0].resize(dim);
    inv_mat[1].resize(dim);
    T deter = recursive_det(B,dim);
    inv_mat[0][0] = B[1][1]/deter;
    inv_mat[0][1] = -B[0][1]/deter;
    inv_mat[1][0] = -B[1][0]/deter;
    inv_mat[1][1] = B[0][0]/deter;
  }else {
    T deter = recursive_det(B,dim);
    inv_mat.resize(dim);
    for(int i=0; i<dim; i++){
      inv_mat[i].resize(dim);
      for(int j=0; j<dim; j++){
        inv_mat[i][j] = pow(-1.0,(i+j)%2)*recursive_det(get_cofactor(B,dim,j,i),dim-1)/deter;
      }
    }
  }
  return inv_mat; 
}

template<class T>
vector<vector<T> > matrix<T>::Leverrier_inv(){
  vector<vector<T> > in;
  if(!is_square()){
    cout << "The matrix is rectangular. This method is valid only for square matrices." << endl;
  }else{
    vector<vector<T> > B = A;
    identity();
    vector<vector<T> > I = A;
    A = B;
    B = I;
    T alpha;
    for(unsigned int i=0; i<N; i++){
      alpha = trace(prod(A,B))/(i+1);
      if(i<N-1){
        B = prod(A,B);
        for(unsigned int j=0; j<N; j++){
          for(unsigned int k=0; k<M; k++){
            if(j==k){
              B[j][k] = -B[j][k] + alpha;
            }else{
              B[j][k] = -B[j][k];
            }
          }
        }
      }
    }
    if(abs(alpha)==0){
      cout << "Alpha at step N is equal to zero. Inverse cannot be computed and left to zero." << endl;
    }else{
      in.resize(N);
      for(unsigned int i=0; i<N; i++){
        in[i].resize(M);
        for(unsigned int j=0; j<M; j++){
          in[i][j] = B[i][j]/alpha;
        }
      }
      in = B/alpha;
    }
  }
  return in;
}

template<class T>
vector<vector<T> > matrix<T>::LU_inv(){
  LU_pp();
  vector<vector<T> > in;
  in.resize(N);
  for(unsigned int i=0; i<N; i++){
    in.resize(M);
  }
  vector<T> x, y, b, ei;
  ei.resize(N);
  for(unsigned int i=0; i<N; i++){
    ei[i] = 0;
  }
  for(unsigned int i=0; i<M; i++){
    if(i>0){
      ei[i-1] = 0;
      ei[i]  = 1;
    }else{
      ei[i]  = 1;
    }
    b.clear();
    b.resize(N);
    for(unsigned int j=0; j<N; j++){
      b[j] = 0;
      for(unsigned int k=0; k<M; k++){
        b[j] += P[j][k]*ei[k];
      }
    }
    y = forw_subs(L,b);
    x = back_subs(U,y);
    for(unsigned int k=0; k<N; k++){
      in[k][i] = x[k];
    }
  }
  return in;
}

template<class T>
void matrix<T>::det(int flag){
  if(!is_square()){
    cout << "The matrix is rectangular. This method is valid only for square matrices." << endl;
    exit(0);
  }else{
    switch(flag){
      case 0:
        d = LU_det();
        break;
      case 1:
        d = recursive_det(A,N);
        break;
    }
  }
}

template<class T>
void matrix<T>::inv(int flag){
  det();
  if(abs(d)==0){
    cout << "The matrix is singular." << endl;
    exit(0);
  }else{
    switch(flag){
      case 0:
        inverse = LU_inv();
        break;
      case 1:
        inverse = Leverrier_inv();
        break;
      case 2:
        inverse = recursive_inv(A,N);
        break;
    }
  }
}

template<class T>
bool matrix<T>::is_zero(){
  for(unsigned int i=0; i<N; i++){
    for(unsigned int j=0; j<M; j++){
      if(A[i][j]!=0){
        return false;
      }
    }
  }
  return true;
}

template<class T>
bool matrix<T>::is_zero(matrix<T> mat){
  vector<int> sizes = mat.size();
  for(unsigned int i=0; i<sizes[0]; i++){
    for(unsigned int j=0; j<sizes[1]; j++){
      if(mat(i,j)!=0){
        return false;
      }
    }
  }
  return true;
}

template<class T>
bool matrix<T>::is_zero(vector<vector<T> > mat){
  int n = mat.size();
  int m = mat[0].size();
  for(unsigned int i=0; i<n; i++){
    for(unsigned int j=0; j<m; j++){
      if(mat[i][j]!=0){
        return false;
      }
    }
  }
  return true;
}

template<class T>
bool matrix<T>::is_square(){
  if(A.size()==A[0].size()){
    return true;
  }else{
    return false;
  }
}                                                  

template<class T>
bool matrix<T>::is_square(matrix<T> mat){
  vector<int> sizes = mat.size();
  if(sizes[0]==sizes[1]){
    return true;
  }else{
    return false;
  }
}

template<class T>
bool matrix<T>::is_square(vector<vector<T> > mat){
  if(mat.size()==mat[0].size()){
    return true;
  }else{
    return false;
  }
}

template<class T>
bool matrix<T>::is_symm(){
  if(!is_square()){
    cout << "The matrix is rectangular." << endl;
    return false;
  }
  for(unsigned int i=0; i<N; i++){
    for(unsigned int j=0; j<M; j++){
      if(j<i){
        if(A[i][j]!=A[j][i]){
          return false;
        }
      }
    }
  }
  return true;
}         

template<class T>
int matrix<T>::is_posdef(){
  if(!is_square()){
    cout << "The matrix is rectangular." << endl;
    return false;
  }else{
    return true;
  }
}

//-----------------
// API functions 
//-----------------

template<class T>
vector<int> matrix<T>::size(){
  vector<int> result;
  result.resize(2);
  result[0] = N;
  result[1] = M;
  return result;
}

template<class T>
T matrix<T>::provide_det(){
  return d;
}

template<class T>
vector<vector<T> > matrix<T>::provide_inv(){
  return inverse;
}

template<class T>
vector<vector<T> > matrix<T>::provide_hermitian(){
  return H;
}

template<class T>
vector<vector<T> > matrix<T>::provide_L(){
  return L;
}

template<class T>
vector<vector<T> > matrix<T>::provide_U(){
  return U;
}

template<class T>
vector<vector<T> > matrix<T>::provide_P(){
  return P;
}

template<class T>
unsigned int matrix<T>::provide_floatop(){
  return floats;
}

template<class T>
string matrix<T>::matlab_writer(string name){
  ostringstream output;
  output << name << " = [";
  if(is_same<T,complex<double> >::value){
    for(int i=0; i<N; i++){
      if(i!=0){
        output << "    ";
      }
      for(int j=0; j<M; j++){
        if(imag(A[i][j])<0.0){
          output << real(A[i][j]) << imag(A[i][j]) << "i";
        }else{
          output << real(A[i][j]) << "+" << imag(A[i][j]) << "i";
        }
        if(j==M-1 && i!=N-1){
          output << "; ..." << endl;
        }else if (j==M-1 && i==N-1){
          output << "];" << endl;
        }else{
          output << ", ";
        }
      }
    }
  }else{
    for(int i=0; i<N; i++){
      if(i!=0){
        output << "    ";
      }
      for(int j=0; j<M; j++){
        output << A[i][j];
        if(j==M-1 && i!=N-1){
          output << "; ..." << endl;
        }else if (j==M-1 && i==N-1){
          output << "];" << endl;
        }else{
          output << ", ";
        }
      }
    }
  }
  return output.str();
}

//--------------------------
// Overloading of operators 
//--------------------------

template<class T>
void matrix<T>::operator=(const matrix<T>& rhs){
  N = rhs.N;
  M = rhs.M;
  A.resize(N);
  for(unsigned int i=0; i<N; i++){
    A[i].resize(M);
    for(unsigned int j=0; j<N; j++){
      A[i][j] = rhs.A[i][j];
    }
  }
}

template<class U>
matrix<U> operator+(const matrix<U>& lhs, const matrix<U>& rhs){
  if(lhs.N!=rhs.N || lhs.M!=rhs.M){
    cout << "Error: dimensions not matching" << endl;
    exit(0);
  }
  matrix<U> mat(lhs.N,lhs.M);
  for(unsigned int i=0; i<lhs.N; i++){
    for(unsigned int j=0; j<lhs.M; j++){
      mat.A[i][j] = lhs.A[i][j] + rhs.A[i][j];
    }
  }
  return mat;
}

template<class U>
matrix<U> operator-(const matrix<U>& lhs, const matrix<U>& rhs){
  if(lhs.N!=rhs.N || lhs.M!=rhs.M){
    cout << "Error: dimensions not matching" << endl;
    exit(0);
  }
  matrix<U> mat(lhs.N,lhs.M);
  for(unsigned int i=0; i<lhs.N; i++){
    for(unsigned int j=0; j<lhs.M; j++){
      mat.A[i][j] = lhs.A[i][j] - rhs.A[i][j];
    }
  }
  return mat;
}

template<class U>
matrix<U> operator*(const matrix<U>& lhs, const matrix<U>& rhs){
  if(lhs.N!=rhs.N || lhs.M!=rhs.M){
    cout << "Error: dimensions not matching" << endl;
    exit(0);
  }
  matrix<U> mat(lhs.N,lhs.M);
  for(unsigned int i=0; i<lhs.N; i++){
    for(unsigned int j=0; j<lhs.M; j++){
      mat.A[i][j] = lhs.A[i][j] * rhs.A[i][j];
    }
  }
  return mat;
}

template<class U>
matrix<U> operator/(const matrix<U>& lhs, const matrix<U>& rhs){
  if(lhs.N!=rhs.N || lhs.M!=rhs.M){
    cout << "Error: dimensions not matching" << endl;
    exit(0);
  }
  matrix<U> mat(lhs.N,lhs.M);
  for(unsigned int i=0; i<lhs.N; i++){
    for(unsigned int j=0; j<lhs.M; j++){
      if(abs(rhs.A[i][j])!=0){
        mat.A[i][j] = lhs.A[i][j] / rhs.A[i][j];
      }else{
        cout << "Error: element (" << i << ", " << j << ") of right hand side matrix is equal to zero. Result is left to zero. Check the computation."  << endl;
      }
    }
  }
  return mat;
}

template<class U>
matrix<U> operator*(const U& lhs, const matrix<U>& rhs){
  matrix<U> mat(rhs.N,rhs.M);
  for(unsigned int i=0; i<rhs.N; i++){
    for(unsigned int j=0; j<rhs.M; j++){
      mat.A[i][j] = rhs.A[i][j] * lhs;
    }
  }
  return mat;
}

template<class U>
matrix<U> operator/(const matrix<U>& lhs, const U& rhs){
  matrix<U> mat(lhs.N,lhs.M);
  for(unsigned int i=0; i<lhs.N; i++){
    for(unsigned int j=0; j<lhs.M; j++){
      mat.A[i][j] = lhs.A[i][j] / rhs;
    }
  }
  return mat;
}

template<class T>
void matrix<T>::operator+=(const matrix<T>& rhs){
  if(N!=rhs.N || M!=rhs.M){
    cout << "Error: dimensions not matching" << endl;
    exit(0);
  }
  for(unsigned int i=0; i<N; i++){
    for(unsigned int j=0; j<M; j++){
      A[i][j] += rhs.A[i][j];
    }
  }
}

template<class T>
void matrix<T>::operator-=(const matrix<T>& rhs){
  if(N!=rhs.N || M!=rhs.M){
    cout << "Error: dimensions not matching" << endl;
    exit(0);
  }
  for(unsigned int i=0; i<N; i++){
    for(unsigned int j=0; j<M; j++){
      A[i][j] -= rhs.A[i][j];
    }
  }
}

template<class T>
void matrix<T>::operator*=(const matrix<T>& rhs){
  if(N!=rhs.N || M!=rhs.M){
    cout << "Error: dimensions not matching" << endl;
    exit(0);
  }
  for(unsigned int i=0; i<N; i++){
    for(unsigned int j=0; j<M; j++){
      A[i][j] *= rhs.A[i][j];
    }
  }
}

template<class T>
void matrix<T>::operator/=(const matrix<T>& rhs){
  if(N!=rhs.N || M!=rhs.M){
    cout << "Error: dimensions not matching" << endl;
    exit(0);
  }
  for(unsigned int i=0; i<N; i++){
    for(unsigned int j=0; j<M; j++){
      A[i][j] /= rhs.A[i][j];
    }
  }
}

template<class T>
void matrix<T>::operator*=(const T& rhs){
  if(N!=rhs.N || M!=rhs.M){
    cout << "Error: dimensions not matching" << endl;
    exit(0);
  }
  for(unsigned int i=0; i<N; i++){
    for(unsigned int j=0; j<M; j++){
      A[i][j] *= rhs;
    }
  }
}

template<class T>
void matrix<T>::operator/=(const T& rhs){
  if(N!=rhs.N || M!=rhs.M){
    cout << "Error: dimensions not matching" << endl;
    exit(0);
  }
  for(unsigned int i=0; i<N; i++){
    for(unsigned int j=0; j<M; j++){
      A[i][j] /= rhs;
    }
  }
}

template<class U>
bool operator==(const matrix<U>& lhs, const matrix<U>& rhs){
  if(lhs.N!=rhs.N || lhs.M!=rhs.M){
    cout << "Error: dimensions not matching" << endl;
    return false;
  }
  for(unsigned int i=0; i<lhs.N; i++){
    for(unsigned int j=0; j<lhs.M; j++){
      if(lhs.A[i][j] != rhs.A[i][j]){
        return false;
      }
    }
  }
  return true;
}

template<class U>
bool operator!=(const matrix<U>& lhs, const matrix<U>& rhs){
  if(lhs.N!=rhs.N || lhs.M!=rhs.M){
    cout << "Error: dimensions not matching" << endl;
    return true;
  }
  for(unsigned int i=0; i<lhs.N; i++){
    for(unsigned int j=0; j<lhs.M; j++){
      if(lhs.A[i][j] != rhs.A[i][j]){
        return true;
      }
    }
  }
  return false;
}

template<class T>
T &matrix<T>::operator()(int i, int j){
  if(i>N && j>M){
    cout << "Error: row and column indeces out of bound" << endl;
    exit(0);
  }
  else if(i>N){
    cout << "Error: row index out of bound" << endl;
    exit(0);
  }
  else if(j>M){
    cout << "Error: column index out of bound" << endl;
    exit(0);
  }else if(N==1){
    return A[j][i];
  }else{
    return A[i][j];
  }
}

template<class T>
matrix<T> &matrix<T>::operator()(int i1, int i2, int j1, int j2){
  if((i1>N || i2>N) && (j1>M || j2>M)){
    cout << "Error: row and column indeces out of bound" << endl;
    exit(0);
  }
  else if(i1>N || i2>N){
    cout << "Error: row indeces out of bound" << endl;
    exit(0);
  }
  else if(j1>M || j2>M){
    cout << "Error: column indeces out of bound" << endl;
    exit(0);
  }else{
    vector<vector<T> > submatrix;
    submatrix.resize(i2-i1+1);
    for(unsigned int i=i1; i<=i2; i++){
      submatrix[i-i1].resize(i2-i1+1);
      for(unsigned int j=j1; j<=j2; j++){
        submatrix[i-i1][j-j1] = A[i][j];
      }
    }
    matrix<T> res(submatrix);
    return res;
  }
}

template<class U>
ostream& operator<<(ostream &output, const matrix<U>& rhs){
  if(is_same<U,complex<double> >::value){
    for(int i=0; i<rhs.N; i++){
      for(int j=0; j<rhs.M; j++){
        if(imag(rhs.A[i][j])<0.0){
          output << real(rhs.A[i][j]) << imag(rhs.A[i][j]) << "i";
        }else{
          output << real(rhs.A[i][j]) << "+" << imag(rhs.A[i][j]) << "i";
        }
        if(j==rhs.M-1){
          output << endl;
        }else{
          output << ", ";
        }
      }
    }
  }else{
    for(int i=0; i<rhs.N; i++){
      for(int j=0; j<rhs.M; j++){
        output << rhs.A[i][j];
        if(j==rhs.M-1){
          output << endl;
        }else{
          output << ", ";
        }
      }
    }
  }
  return output;
}

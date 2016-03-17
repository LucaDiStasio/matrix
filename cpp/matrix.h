/*
 * =====================================================================================
 * 
 * Copyright (c) 2016 Université de Lorraine & Luleå tekniska universitet
 * Author: Luca Di Stasio <luca.distasio@gmail.com>
 *                        <luca.distasio@ingpec.eu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * =====================================================================================
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <random>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <typeinfo>
#include <type_traits> // C++0x
//#include <tr1/type_traits> // C++03, use std::tr1
#include <vector>

using namespace std;

//============================================================================//
//============================================================================//
/*
           A class to implement matrix operations and matrix algebra
*/
//============================================================================//
//============================================================================//


//===================================================
//==================  HEADER  =======================
//===================================================

template<class T = double>
class matrix {

  //===================================================  
  //                  Variables
  //===================================================
  private:

    // Control parameters
    bool measure;                                    // True --> count floating point and integer operations, false otherwise
    unsigned int floats;                             // Number of floating point operations performed
    unsigned int ips;                                // Number of instructions (integer operations) performed

  // Output quantities
    int N;                                           // Number of rows of A
    int M;                                           // Number of columns of A
    vector<vector<T> > A;                            // Matrix
    vector<vector<T> > L;                            // Lower triangular matrix (for LU decomposition)
    vector<vector<T> > U;                            // Upper triangular matrix (for LU decomposition)
    vector<vector<T> > P;                            // Permutation matrix  (for LU decomposition)
    vector<vector<T> > Q;                            // Permutation matrix  (for LU decomposition with full pivoting)
    T d;                                             // Determinant
    vector<vector<T> > inverse;                      // Inverse matrix
    vector<vector<T> > H;                            // Hermitian of the matrix (transpose and complex conjugate of each element)
    vector<T> invariants;                            // Matrix invariants

  
  //===================================================  
  //                      Methods
  //===================================================  
  public:
  
    // Constructor (default)
    matrix();
    
    // Constructor (number of rows for square matrix)
    matrix(int rows);
    
    // Constructor (number of rows for square matrix, measure performance)
    matrix(int rows, bool meas);
    
    // Constructor (number of rows and columns for rectangular matrix)
    matrix(int rows, int columns);
    
    // Constructor (number of rows and columns for rectangular matrix, measure performance)
    matrix(int rows, int columns, bool meas);
    
    // Constructor (square matrix)
    matrix(T inp[], int rows);
    
    // Constructor (square matrix, measure performances)
    matrix(T inp[], int rows, bool meas);
    
    // Constructor (rectangular matrix)
    matrix(T inp[], int rows, int columns);
    
    // Constructor (rectangular matrix, measure performances)
    matrix(T inp[], int rows, int columns, bool meas);
    
    // Constructor (matrix)
    matrix(vector<vector<T> > inp);
    
    // Constructor (matrix, measure performances)
    matrix(vector<vector<T> > inp, bool meas);
    
    //Destructor
    ~matrix();
  
    // General tools
    void clear();                                                                         // Clear matrix
    
    // Linear algebra
    void init_random(double a=0.0, double b=1.0);                                         // Generate a random matrix
  
    void identity();                                                                      // Set A equal to the identity matrix
    void identity(T value);                                                               // Set A equal to value*identity matrix
    
    void ones();                                                                          // Set all elements of A equal to 1
    void ones(T value);                                                                   // Set all elements of A equal to value
    
    void diag(T vec[], int length, int pos = 0);                                          // Set the element of the pos-th diagonal of A (0 -> main diagonal, +pos -> upper diagonal, -pos -> lower diagonal)
    void diag(vector<T> vec, int pos = 0);                                                // Set the element of the pos-th diagonal of A (0 -> main diagonal, +pos -> upper diagonal, -pos -> lower diagonal)
    
    void hermitian();                                                                     // Compute Hermitian matrix (transpose and complex conjugate)
     
    T trace();                                                                            // Trace of A
    T trace(vector<vector<T> > mat);                                                      // Trace of mat
    T trace(matrix<T> mat);                                                               // Trace of mat
    
    T matnorm1();                                                                         // Matrix 1-norm of A (if A is a matrix)
    T matnormInf();                                                                       // Matrix Infinity-norm of A  (if A is a matrix)
    T matnormFrob();                                                                      // Matrix Frobenius-norm of A  (if A is a matrix)
    
    T vecnorm2();                                                                         // Vector 2-norm of A  (if A is a vector)
    
    vector<vector<T> > absv();                                                            // Absolute value of matrix (matrix of absolute values of elements)
    vector<vector<T> > absv(vector<vector<T> > mat);                                      // Absolute value of matrix (matrix of absolute values of elements)
    vector<vector<T> > absv(matrix<T> mat);                                               // Absolute value of matrix (matrix of absolute values of elements)

    matrix<T> absm();                                                                     // Absolute value of matrix (matrix of absolute values of elements)
    matrix<T> absm(vector<vector<T> > mat);                                               // Absolute value of matrix (matrix of absolute values of elements)
    matrix<T> absm(matrix<T> mat);                                                        // Absolute value of matrix (matrix of absolute values of elements)
    
    vector<int> min();                                                                    // Return indices of minimum value in A
    vector<int> min(vector<vector<T> > mat);                                              // Return indices of minimum value in mat
    vector<int> min(matrix mat);                                                          // Return indices of minimum value in mat
    
    vector<int> max();                                                                    // Return indices of maximum value in A
    vector<int> max(vector<vector<T> > mat);                                              // Return indices of maximum value in mat
    vector<int> max(matrix<T> mat);                                                       // Return indices of maximum value in mat
    
    matrix<T> vec_min(int flag);                                                          // Return minimum values over columns (flag = 0, default) and return a row vector or over rows (flag = 1) and return a column vector of A
    matrix<T> vec_min(vector<vector<T> > mat, int flag);                                  // Return minimum values over columns (flag = 0, default) and return a row vector or over rows (flag = 1) and return a column vector of mat
    matrix<T> vec_min(matrix<T> mat, int flag);                                           // Return minimum values over columns (flag = 0, default) and return a row vector or over rows (flag = 1) and return a column vector of mat
    
    matrix<T> vec_max(int flag);                                                          // Return maximum values over columns (flag = 0, default) and return a row vector or over rows (flag = 1) and return a column vector of A
    matrix<T> vec_max(vector<vector<T> > mat, int flag);                                  // Return maximum values over columns (flag = 0, default) and return a row vector or over rows (flag = 1) and return a column vector of mat
    matrix<T> vec_max(matrix<T> mat, int flag);                                           // Return maximum values over columns (flag = 0, default) and return a row vector or over rows (flag = 1) and return a column vector of mat
    
    matrix<T> right_prod(matrix<T> mat);                                                  // Compute the matrix product A*mat
    matrix<T> left_prod(matrix<T> mat);                                                   // Compute the matrix product mat*A
    matrix<T> prod(matrix<T> mat1, matrix<T> mat2);                                       // Compute the matrix product mat1*mat2
    vector<vector<T> > prod(vector<vector<T> > mat1, vector<vector<T> > mat2);            // Compute the matrix product mat1*mat2
    
    T right_scalar_prod(matrix<T> mat);                                                   // Compute the scalar product between vectors A*mat (represented as N x 1 or 1 x N matrices)
    T left_scalar_prod(matrix<T> mat);                                                    // Compute the scalar product between vectors mat*A (represented as N x 1 or 1 x N matrices)
    T scalar_prod(matrix<T> mat1, matrix<T> mat2);                                        // Compute the scalar product between vectors mat1*mat2 (represented as N x 1 or 1 x N matrices)
    
    matrix<T> right_tensor_prod(matrix<T> mat);                                           // Compute the tensor product A*mat
    matrix<T> left_tensor_prod(matrix<T> mat);                                            // Compute the tensor product mat*A
    matrix<T> tensor_prod(matrix<T> mat1, matrix<T> mat2);                                // Compute the tensor product mat1*mat2
    
    vector<vector<T> > get_cofactor(vector<vector<T> > B, int n, int row, int col);       // Compute the minor
    
    matrix<T> forw_subs(matrix<T> b);                                                     // Forward substitutions using upper-triangular matrix U and input vector
    matrix<T> forw_subs(matrix<T> Linp, matrix<T> b);                                     // Forward substitutions using input upper-triangular matrix and input vector
    vector<T> forw_subs(vector<T> b);                                                     // Forward substitutions using upper-triangular matrix U and input vector
    vector<T> forw_subs(vector<vector<T> > Linp, vector<T> b);                            // Forward substitutions using input upper-triangular matrix and input vector
    
    matrix<T> back_subs(matrix<T> b);                                                     // Backward substitutions using lower-triangular matrix L and input vector
    matrix<T> back_subs(matrix<T> Uinp, matrix<T> b);                                     // Backward substitutions using input lower-triangular matrix and input vectors
    vector<T> back_subs(vector<T> b);                                                     // Backward substitutions using lower-triangular matrix L and input vector
    vector<T> back_subs(vector<vector<T> > Uinp, vector<T> b);                            // Backward substitutions using input lower-triangular matrix and input vectors
    
    void LU();                                                                            // Compute LU factorization without pivoting - Doolittle version, i.e. L_ii = 1
    void LU_pp();                                                                         // Compute LU factorization with partial pivoting
    void LU_fp();                                                                         // Compute LU factorization with full pivoting
    void Cholesky();                                                                      // Compute Cholesky factorization
    
    T recursive_det(vector<vector<T> > B, int n);                                         // Recursive computation of determinant
    T LU_det();                                                                           // Computation of determinant with LU factorization
    
    vector<vector<T> > recursive_inv(vector<vector<T> > B, int dim);                      // Recursive computation of inverse
    vector<vector<T> > Leverrier_inv();                                                   // Computation of inverse with Leverrier method
    vector<vector<T> > LU_inv();                                                          // Computation of inverse with LU factorization
    
    void det(int flag = 0);                                                               // Compute determinant (flag = 0 --> LU, default; flag = 1 --> Leverrier; flag = 2 --> recursive)
    void inv(int flag = 0);                                                               // Compute inverse (flag = 0 --> LU, default; flag = 1 --> Leverrier; flag = 2 --> recursive)
   
    bool is_zero();                                                                       // True if all elements of A are equal to zero
    bool is_zero(matrix<T> mat);                                                          // True if all elements of mat are equal to zero
    bool is_zero(vector<vector<T> > mat);                                                 // True if all elements of mat are equal to zero
    
    bool is_square();                                                                     // True if A is a square matrix
    bool is_square(matrix<T> mat);                                                        // True if mat is a square matrix
    bool is_square(vector<vector<T> > mat);                                               // True if mat is a square matrix
    
    bool is_symm();                                                                       // True if A is symmetric, false otherwise
    int is_posdef();                                                                      // 1 if A is positive definite, -1 if it's negative definite, 0 if it's indefinite
    
    vector<int> size();                                                                   // Provide size of matrix [N x M] = [rows x cols]

    T provide_det();                                                                      // Provide determinant to external program
    vector<vector<T> > provide_inv();                                                     // Provide inverse matrix to external program
    vector<vector<T> > provide_hermitian();                                               // Provide transpose matrix to external program
    vector<vector<T> > provide_L();                                                       // Provide L matrix (from LU factorization) to external program
    vector<vector<T> > provide_U();                                                       // Provide U matrix (from LU factorization) to external program
    vector<vector<T> > provide_P();                                                       // Provide P matrix (from LU factorization) to external program
    unsigned int provide_floatop();                                                       // Provide total number of floating point operations to external program
    
    string matlab_writer(string name = "A");                                              // Print A in matlab format
    
    // Overloading of operators
    void operator=(const matrix<T> &rhs);                                                 // Set A equal to an external matrix
    
    template<class U>
    friend matrix<U> operator+(const matrix<U>& lhs, const matrix<U>& rhs);               // Sum two matrices element-by-element
    
    template<class U>
    friend matrix<U> operator-(const matrix<U>& lhs, const matrix<U>& rhs);               // Subtract two matrices element-by-element
    
    template<class U>
    friend matrix<U> operator*(const matrix<U>& lhs, const matrix<U>& rhs);               // Multiply two matrices element-by-element
    
    template<class U>
    friend matrix<U> operator/(const matrix<U>& lhs, const matrix<U>& rhs);               // Divide two matrices element-by-element
    
    template<class U>
    friend matrix<U> operator*(const U& lhs, const matrix<U>& rhs);                       // Multiply by scalar element-by-element
    
    template<class U>
    friend matrix<U> operator/(const matrix<U>& lhs, const U& rhs);                       // Divide by scalar element-by-element
    
    void operator+=(const matrix<T>& rhs);                                                // Sum two matrices element-by-element
    void operator-=(const matrix<T>& rhs);                                                // Subtract two matrices element-by-element
    void operator*=(const matrix<T>& rhs);                                                // Multiply two matrices element-by-element
    void operator/=(const matrix<T>& rhs);                                                // Divide two matrices element-by-element
    void operator*=(const T& rhs);                                                        // Multiply by scalar element-by-element
    void operator/=(const T& rhs);                                                        // Divide by scalar element-by-element
    
    template<class U>
    friend bool operator==(const matrix<U>& lhs, const matrix<U>& rhs);                   // Compare two matrices (equality)
    
    template<class U>
    friend bool operator!=(const matrix<U>& lhs, const matrix<U>& rhs);                   // Compare two matrices (inequality)
    
    T &operator()(int i, int j = 0);                                                      // Access element (i,j)
    matrix<T> &operator()(int i1, int i2, int j1, int j2);                                // Access elements (i1:i2,j1:j2)
    
    template<class U>
    friend ostream& operator<<(ostream& output, const matrix<U>& rhs);                    // Print matrix to screen
    
    ostream& matlab_writer(ostream& output, string name);                                 // Print A in matlab format
};


#endif
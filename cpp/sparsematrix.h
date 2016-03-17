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

#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

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
    A class to implement sparse matrix operations and sparse matrix algebra
*/
//============================================================================//
//============================================================================//


//===================================================
//==================  HEADER  =======================
//===================================================

template<class T = double>
class sparsematrix {

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
    int NNZ;                                         // Number of non-zero elements of A
    vector<T> values;                                // Non-zero matrix values
    vector<int> indeces;                             // Indeces of non-zero elements (helical numbering)
    double ratio;                                    // Ratio of number of non-zero elements to total number of elements
    double eta;                                      // Effectiveness eta = 2*NNZ/(N*M), it measures the effectiveness of using 
                                                   // the sparse vs the full matrix format in terms of number of values stored
  
    //===================================================  
    //                      Methods
    //===================================================  
    public:
  
    // Constructor (default)
      sparsematrix();
      
      // Constructor (number of rows for square matrix)
      sparsematrix(int rows);
      
      // Constructor (number of rows for square matrix, measure performance)
      sparsematrix(int rows, bool meas);
      
      // Constructor (number of rows and columns for rectangular matrix)
      sparsematrix(int rows, int columns);
      
      // Constructor (number of rows and columns for rectangular matrix, measure performance)
      sparsematrix(int rows, int columns, bool meas);
      
      // Constructor (square matrix)
      sparsematrix(T inp[], int ind[], int rows);
      
      // Constructor (square matrix, measure performances)
      sparsematrix(T inp[], int ind[], int rows, bool meas);
      
      // Constructor (rectangular matrix)
      sparsematrix(T inp[], int ind[], int rows, int columns);
      
      // Constructor (rectangular matrix, measure performances)
      sparsematrix(T inp[], int ind[], int rows, int columns, bool meas);
      
      // Constructor (Square matrix)
      sparsematrix(vector<T> inp, vector<vector<int> > ind, int rows);
      
      // Constructor (Square matrix, measure performances)
      sparsematrix(vector<T> inp, vector<vector<int> > ind, int rows, bool meas);
      
      // Constructor (Rectangular matrix)
      sparsematrix(vector<T> inp, vector<vector<int> > ind, int rows, int columns);
      
      // Constructor (Rectangular matrix, measure performances)
      sparsematrix(vector<T> inp, vector<vector<int> > ind, int rows, int columns, bool meas);
      
      //Destructor
      ~sparsematrix();
      
      // General tools
      void clear();                                                                         // Clear matrix
      
      // Linear algebra
      void init_random(double a=0.0, double b=1.0);                                         // Generate a random matrix
    
      void identity();                                                                      // Set A equal to the identity matrix
      void identity(T value);                                                               // Set A equal to value*identity matrix
      
      void diag(T vec[], int length, int pos = 0);                                          // Set the element of the pos-th diagonal of A (0 -> main diagonal, +pos -> upper diagonal, -pos -> lower diagonal)
      void diag(vector<T> vec, int pos = 0);                                                // Set the element of the pos-th diagonal of A (0 -> main diagonal, +pos -> upper diagonal, -pos -> lower diagonal)
      
      matrix<T> right_prod(matrix<T> mat);                                                  // Compute the matrix product A*mat
      matrix<T> left_prod(matrix<T> mat);                                                   // Compute the matrix product mat*A
      matrix<T> prod(matrix<T> mat1, matrix<T> mat2);                                       // Compute the matrix product mat1*mat2
      T right_scalar_prod(matrix<T> mat);                                                   // Compute the scalar product between vectors A*mat (represented as N x 1 or 1 x N matrices)
      T left_scalar_prod(matrix<T> mat);                                                    // Compute the scalar product between vectors mat*A (represented as N x 1 or 1 x N matrices)
      T scalar_prod(matrix<T> mat1, matrix<T> mat2);                                        // Compute the scalar product between vectors mat1*mat2 (represented as N x 1 or 1 x N matrices)
      matrix<T> right_tensor_prod(matrix<T> mat);                                           // Compute the tensor product A*mat
      matrix<T> left_tensor_prod(matrix<T> mat);                                            // Compute the tensor product mat*A
      matrix<T> tensor_prod(matrix<T> mat1, matrix<T> mat2);                                // Compute the tensor product mat1*mat2
      
      bool is_zero();                                                                       // True if all elements of A are equal to zero
      bool is_zero(matrix<T> mat);                                                          // True if all elements of mat are equal to zero
      bool is_zero(vector<vector<T> > mat);                                                 // True if all elements of mat are equal to zero
      
      bool is_square();                                                                     // True if A is a square matrix
      bool is_square(matrix<T> mat);                                                        // True if mat is a square matrix
      bool is_square(vector<vector<T> > mat);                                               // True if mat is a square matrix
      
      bool is_symm();                                                                       // True if A is symmetric, false otherwise
      int is_posdef();                                                                      // 1 if A is positive definite, -1 if it's negative definite, 0 if it's indefinite
      
      vector<int> size();                                                                   // Provide size of matrix [N x M] = [rows x cols]
      
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
      
      template<class U>
      friend ostream& operator<<(ostream& output, const matrix<U>& rhs);                    // Print matrix to screen
      
      ostream& matlab_writer(ostream& output, string name);                                 // Print A in matlab format
};


#endif
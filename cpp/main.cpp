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

//#include <boost/program_options/options_description.hpp>
//#include <boost/program_options/parsers.hpp>
//#include <boost/program_options/variables_map.hpp>

#include <iostream>
#include <fstream>

#include "datetime.h"
#include "matrix.h"
#include "matrix.cpp"

using namespace std;
//using namespace boost;
//using namespace boost::program_options;

//==================  MAIN  =======================
int main(int argc, char** argv){
    
    string WorkDir = "/home/ubuntu/workspace/Matrix";
    string OutDir = "Code";
    
    bool test1 = false;
    bool test2 = true;
    
    bool print_screen = true;
    bool print_dat = true;
    bool print_csv = true;
    bool print_mat = true;
    
    datetime today;
    string format = "B M a, Y - H:i:s";
    string date = today.get_customtimestamp(format);
    
    const char* alphabet[] = {"A", "B", "C", "D", "E", "F", "G", "H", "K", "J", "I", "L", "M", "N", "O",  "P",  "Q",  "R",  "S",  "T",  "U",  "V",  "X",  "Y", "W", "Z"};
    vector<string> names(alphabet, alphabet + 26);
  
    if(test1){
        double matrix1[] = {1.04, 57.9, -36.4, 15.5};
        double matrix2[] = {12.69, 9.75, -95.75, 91.33, 27.84, 96.48, 63.23, 54.68, 15.76};
        double matrix3[] = {97.05, -14.18, 95.94, 93.39, 95.71, -42.17, 65.57, 67.87, 48.53, 91.57, 3.57, 75.77, 80.02, 79.22, 84.91, 74.31};
      
        int matrix4[] = {5,8,6,-4};
        int matrix5[] = {56,1,-7,54,18,-61,21,3,9};
        int matrix6[] = {24,-28,19,97,-42,-66,16,84,39,71,-75,57,-33,13,6,47};
      
        complex<double> matrix7[] = {{79.64,65.12},{3.74,34.94},{-93.14,9.07},{74.54,-22.31}};
        complex<double> matrix8[] = {{19.34,-31.97},{36.01,5.81},{-7.18,-30.81},{1.09,-11.47},{12.88,27.07},{8.13,4.15},{-35.48,47.66},{28.23,-91.99},{82.92,13.61}};
        complex<double> matrix9[] = {{11.48,-43.44},{34.17,78.16},{18.73,17.63},{-69.11,48.73},{12.10,-17.9},{-31.17,71.93},{66.54,4.24},{62.02,67.00},{-10.34,-21.41},{-28.57,-9.51},{65.42,25.83},{4.79,14.97},{3.27,-71.4},{48.81,42.79},{33.35,3.80},{39.27,59.44}};
      
        double f1 = 3.5;
        int f2 = 3;
        complex<double> f3(3.5,0.0);
      
        matrix<double> mat1(matrix1,2);
        matrix<double> mat2(matrix2,3);
        matrix<double> mat3(matrix3,4);
        matrix<int> mat4(matrix4,2);
        matrix<int> mat5(matrix5,3);
        matrix<int> mat6(matrix6,4);
        matrix<complex<double> > mat7(matrix7,2);
        matrix<complex<double> > mat8(matrix8,3);
        matrix<complex<double> > mat9(matrix9,4);
    
        mat1.det();
        mat1.inv();
        mat1.hermitian();
        double deter1 = mat1.provide_det();
        matrix<double> inverse1(mat1.provide_inv());
        matrix<double> transpose1(mat1.provide_hermitian());
      
        mat2.det();
        mat2.inv();
        mat2.hermitian();
        double deter2 = mat2.provide_det();
        matrix<double> inverse2(mat2.provide_inv());
        matrix<double> transpose2(mat2.provide_hermitian());
      
        mat3.det();
        mat3.inv();
        mat3.hermitian();
        double deter3 = mat3.provide_det();
        matrix<double> inverse3(mat3.provide_inv());
        matrix<double> transpose3(mat3.provide_hermitian());
      
        mat4.det();
        mat4.inv();
        mat4.hermitian();
        int deter4 = mat4.provide_det();
        matrix<int> inverse4(mat4.provide_inv());
        matrix<int> transpose4(mat4.provide_hermitian());
      
        mat5.det();
        mat5.inv();
        mat5.hermitian();
        int deter5 = mat5.provide_det();
        matrix<int> inverse5(mat5.provide_inv());
        matrix<int> transpose5(mat5.provide_hermitian());
      
        mat6.det();
        mat6.inv();
        mat6.hermitian();
        int deter6 = mat6.provide_det();
        matrix<int> inverse6(mat6.provide_inv());
        matrix<int> transpose6(mat6.provide_hermitian());
      
        mat7.det();
        mat7.inv();
        mat7.hermitian();
        complex<double> deter7 = mat7.provide_det();
        matrix<complex<double> > inverse7(mat7.provide_inv());
        matrix<complex<double> > transpose7(mat7.provide_hermitian());
      
        mat8.det();
        mat8.inv();
        mat8.hermitian();
        complex<double> deter8 = mat8.provide_det();
        matrix<complex<double> > inverse8(mat8.provide_inv());
        matrix<complex<double> > transpose8(mat8.provide_hermitian());
      
        mat9.det();
        mat9.inv();
        mat9.hermitian();
        complex<double> deter9 = mat9.provide_det();
        matrix<complex<double> > inverse9(mat9.provide_inv());
        matrix<complex<double> > transpose9(mat9.provide_hermitian());
        
        vector<matrix<double> > doublemat;
        vector<matrix<int> > intmat;
        vector<matrix<complex<double> > > complexmat;
        doublemat.resize(3);
        intmat.resize(3);
        complexmat.resize(3);
        doublemat[0] = mat1;
        doublemat[1] = mat2;
        doublemat[2] = mat3;
        intmat[0] = mat4;
        intmat[1] = mat5;
        intmat[2] = mat6;
        complexmat[0] = mat7;
        complexmat[1] = mat8;
        complexmat[2] = mat9;
        
        if (print_screen){
            cout << "A = " << endl;
            cout << mat1;
            cout << endl;
            cout << "det(A) = " << deter1 << endl;
            cout << endl;
            cout << "inv(A) = " << endl;
            cout << inverse1;
            cout << endl;
            cout << "A^H = " << endl;
            cout << transpose1;
            cout << endl;
            cout << "3.5*A = " << endl;
            cout << f1*mat1;
            cout << endl;
            cout << "A + inv(A) = " << endl;
            cout << mat1 - inverse1;
            cout << endl;
            cout << "A - A^H = " << endl;
            cout << mat1 -  transpose1;
            cout << endl;
            cout << "A * A^H = " << endl;
            cout << mat1 *  transpose1;
            cout << endl;
            cout << "A / inv(A) = " << endl;
            cout << mat1 /  inverse1;
            cout << endl;
            cout << "A(1,2) = " << endl;
            cout << mat1(0,1) << endl;
            cout << endl;
      
            cout << "B = " << endl;
            cout << mat2;
            cout << endl;
            cout << "det(B) = " << deter2 << endl;
            cout << endl;
            cout << "inv(B) = " << endl;
            cout << inverse2;
            cout << endl;
            cout << "B^H = " << endl;
            cout << transpose2;
            cout << endl;
            cout << "3.5*B = " << endl;
            cout << f1*mat2;
            cout << endl;
            cout << "B + inv(B) = " << endl;
            cout << mat2 - inverse2;
            cout << endl;
            cout << "B - B^H = " << endl;
            cout << mat2 -  transpose2;
            cout << endl;
            cout << "B * B^H = " << endl;
            cout << mat2 *  transpose2;
            cout << endl;
            cout << "B / inv(B) = " << endl;
            cout << mat2 /  inverse2;
            cout << endl;
            cout << "B(1,2) = " << endl;
            cout << mat2(0,1) << endl;
            cout << endl;
      
            cout << "C = " << endl;
            cout << mat3;
            cout << endl;
            cout << "det(C) = " << deter3 << endl;
            cout << endl;
            cout << "inv(C) = " << endl;
            cout << inverse3;
            cout << endl;
            cout << "C^H = " << endl;
            cout << transpose3;
            cout << endl;
            cout << "3.5*C = " << endl;
            cout << f1*mat3;
            cout << endl;
            cout << "C + inv(C) = " << endl;
            cout << mat3 - inverse3;
            cout << endl;
            cout << "C - C^H = " << endl;
            cout << mat3 -  transpose3;
            cout << endl;
            cout << "C * C^H = " << endl;
            cout << mat3 *  transpose3;
            cout << endl;
            cout << "C / inv(C) = " << endl;
            cout << mat3 /  inverse3;
            cout << endl;
            cout << "C(1,2) = " << endl;
            cout << mat3(0,1) << endl;
            cout << endl;
      
            cout << "D = " << endl;
            cout << mat4;
            cout << endl;
            cout << "det(D) = " << deter4 << endl;
            cout << endl;
            cout << "inv(D) = " << endl;
            cout << inverse4;
            cout << endl;
            cout << "D^H = " << endl;
            cout << transpose4;
            cout << endl;
            cout << "3*D = " << endl;
            cout << f2*mat4;
            cout << endl;
            cout << "D + inv(D) = " << endl;
            cout << mat4 - inverse4;
            cout << endl;
            cout << "D - D^H = " << endl;
            cout << mat4 -  transpose4;
            cout << endl;
            cout << "D * D^H = " << endl;
            cout << mat4 *  transpose4;
            cout << endl;
            /*cout << "D / inv(D) = " << endl;
            cout << mat4 /  inverse4;
            cout << endl;*/
            cout << "D(1,2) = " << endl;
            cout << mat4(0,1) << endl;
            cout << endl;
      
            cout << "E = " << endl;
            cout << mat5;
            cout << endl;
            cout << "det(E) = " << deter5 << endl;
            cout << endl;
            cout << "inv(E) = " << endl;
            cout << inverse5;
            cout << endl;
            cout << "E^H = " << endl;
            cout << transpose5;
            cout << endl;
            cout << "3*E = " << endl;
            cout << f2*mat5;
            cout << endl;
            cout << "E + inv(E) = " << endl;
            cout << mat5 - inverse5;
            cout << endl;
            cout << "E - E^H = " << endl;
            cout << mat5 -  transpose5;
            cout << endl;
            cout << "E * E^H = " << endl;
            cout << mat5 *  transpose5;
            cout << endl;
            /*cout << "E / inv(E) = " << endl;
            cout << mat5 /  inverse5;
            cout << endl;*/
            cout << "E(1,2) = " << endl;
            cout << mat5(0,1) << endl;
            cout << endl;
      
            cout << "F = " << endl;
            cout << mat6;
            cout << endl;
            cout << "det(F) = " << deter6 << endl;
            cout << endl;
            cout << "inv(F) = " << endl;
            cout << inverse6;
            cout << endl;
            cout << "F^H = " << endl;
            cout << transpose6;
            cout << endl;
            cout << "3*F = " << endl;
            cout << f2*mat6;
            cout << endl;
            cout << "F + inv(F) = " << endl;
            cout << mat6 - inverse6;
            cout << endl;
            cout << "F - F^H = " << endl;
            cout << mat6 -  transpose6;
            cout << endl;
            cout << "F * F^H = " << endl;
            cout << mat6 *  transpose6;
            cout << endl;
            /*cout << "F / inv(F) = " << endl;
            cout << mat6 /  inverse6;
            cout << endl;*/
            cout << "F(1,2) = " << endl;
            cout << mat6(0,1) << endl;
            cout << endl;
      
            cout << "G = " << endl;
            cout << mat7;
            cout << endl;
            cout << "det(G) = " << real(deter7) << "+" << imag(deter7) << "i" << endl;
            cout << endl;
            cout << "inv(G) = " << endl;
            cout << inverse7;
            cout << endl;
            cout << "G^H = " << endl;
            cout << transpose7;
            cout << endl;
            cout << "3.5*G = " << endl;
            cout << f3*mat7;
            cout << endl;
            cout << "G + inv(G) = " << endl;
            cout << mat7 - inverse7;
            cout << endl;
            cout << "G - G^H = " << endl;
            cout << mat7 -  transpose7;
            cout << endl;
            cout << "G * G^H = " << endl;
            cout << mat7 *  transpose7;
            cout << endl;
            cout << "G / inv(G) = " << endl;
            cout << mat7 /  inverse7;
            cout << endl;
            cout << "G(1,2) = " << endl;
            cout << real(mat7(0,1)) << "+" << imag(mat7(0,1)) << "i" << endl;
            cout << endl;
      
            cout << "H = " << endl;
            cout << mat8;
            cout << endl;
            cout << "det(H) = " <<  real(deter8) << "+" << imag(deter8) << "i" << endl;
            cout << endl;
            cout << "inv(H) = " << endl;
            cout << inverse8;
            cout << endl;
            cout << "H^H = " << endl;
            cout << transpose8;
            cout << endl;
            cout << "3.5*H = " << endl;
            cout << f3*mat8;
            cout << endl;
            cout << "H + inv(H) = " << endl;
            cout << mat8 - inverse8;
            cout << endl;
            cout << "H - H^H = " << endl;
            cout << mat8 -  transpose8;
            cout << endl;
            cout << "H * H^H = " << endl;
            cout << mat8 *  transpose8;
            cout << endl;
            cout << "H / inv(H) = " << endl;
            cout << mat8 /  inverse8;
            cout << endl;
            cout << "H(1,2) = " << endl;
            cout << real(mat8(0,1)) << "+" << imag(mat8(0,1)) << "i" << endl;
            cout << endl;
      
            cout << "I = " << endl;
            cout << mat9;
            cout << endl;
            cout << "det(I) = " <<  real(deter9) << "+" << imag(deter9) << "i" << endl;
            cout << endl;
            cout << "inv(I) = " << endl;
            cout << inverse9;
            cout << endl;
            cout << "I^H = " << endl;
            cout << transpose9;
            cout << endl;
            cout << "3.5*I = " << endl;
            cout << f3*mat9;
            cout << endl;
            cout << "I + inv(I) = " << endl;
            cout << mat9 - inverse9;
            cout << endl;
            cout << "I - I^H = " << endl;
            cout << mat9 -  transpose9;
            cout << endl;
            cout << "I * I^H = " << endl;
            cout << mat9 *  transpose9;
            cout << endl;
            cout << "I / inv(I) = " << endl;
            cout << mat9 /  inverse9;
            cout << endl;
            cout << "I(1,2) = " << endl;
            cout << real(mat9(0,1)) << "+" << imag(mat9(0,1)) << "i" << endl;
            cout << endl;
      }
        
        if (print_dat){
            string OutFileName = "MatrixTest1.dat";
            string OutFilePathString = WorkDir + "/" + OutDir + "/" + OutFileName;
            const char * OutFilePath = OutFilePathString.c_str();
            ofstream outfile;
            outfile.open(OutFilePath);
            if(outfile.is_open()){
                outfile << "%% A = " << endl;
                outfile << mat1;
                outfile << endl;
                outfile << "%% det(A) = " << deter1 << endl;
                outfile << endl;
                outfile << "%% inv(A) = " << endl;
                outfile << inverse1;
                outfile << endl;
                outfile << "%% A^H = " << endl;
                outfile << transpose1;
                outfile << endl;
                outfile << "%% 3.5*A = " << endl;
                outfile << f1*mat1;
                outfile << endl;
                outfile << "%% A + inv(A) = " << endl;
                outfile << mat1 - inverse1;
                outfile << endl;
                outfile << "%% A - A^H = " << endl;
                outfile << mat1 -  transpose1;
                outfile << endl;
                outfile << "%% A * A^H = " << endl;
                outfile << mat1 *  transpose1;
                outfile << endl;
                outfile << "%% A / inv(A) = " << endl;
                outfile << mat1 /  inverse1;
                outfile << endl;
                outfile << "%% A(1,2) = " << endl;
                outfile << mat1(0,1) << endl;
                outfile << endl;
      
                outfile << "%% B = " << endl;
                outfile << mat2;
                outfile << endl;
                outfile << "%% det(B) = " << deter2 << endl;
                outfile << endl;
                outfile << "%% inv(B) = " << endl;
                outfile << inverse2;
                outfile << endl;
                outfile << "%% B^H = " << endl;
                outfile << transpose2;
                outfile << endl;
                outfile << "%% 3.5*B = " << endl;
                outfile << f1*mat2;
                outfile << endl;
                outfile << "%% B + inv(B) = " << endl;
                outfile << mat2 - inverse2;
                outfile << endl;
                outfile << "%% B - B^H = " << endl;
                outfile << mat2 -  transpose2;
                outfile << endl;
                outfile << "%% B * B^H = " << endl;
                outfile << mat2 *  transpose2;
                outfile << endl;
                outfile << "%% B / inv(B) = " << endl;
                outfile << mat2 /  inverse2;
                outfile << endl;
                outfile << "%% B(1,2) = " << endl;
                outfile << mat2(0,1) << endl;
                outfile << endl;
      
                outfile << "%% C = " << endl;
                outfile << mat3;
                outfile << endl;
                outfile << "%% det(C) = " << deter3 << endl;
                outfile << endl;
                outfile << "%% inv(C) = " << endl;
                outfile << inverse3;
                outfile << endl;
                outfile << "%% C^H = " << endl;
                outfile << transpose3;
                outfile << endl;
                outfile << "%% 3.5*C = " << endl;
                outfile << f1*mat3;
                outfile << endl;
                outfile << "%% C + inv(C) = " << endl;
                outfile << mat3 - inverse3;
                outfile << endl;
                outfile << "%% C - C^H = " << endl;
                outfile << mat3 -  transpose3;
                outfile << endl;
                outfile << "%% C * C^H = " << endl;
                outfile << mat3 *  transpose3;
                outfile << endl;
                outfile << "%% C / inv(C) = " << endl;
                outfile << mat3 /  inverse3;
                outfile << endl;
                outfile << "%% C(1,2) = " << endl;
                outfile << mat3(0,1) << endl;
                outfile << endl;
      
                outfile << "%% D = " << endl;
                outfile << mat4;
                outfile << endl;
                outfile << "%% det(D) = " << deter4 << endl;
                outfile << endl;
                outfile << "%% inv(D) = " << endl;
                outfile << inverse4;
                outfile << endl;
                outfile << "%% D^H = " << endl;
                outfile << transpose4;
                outfile << endl;
                outfile << "%% 3*D = " << endl;
                outfile << f2*mat4;
                outfile << endl;
                outfile << "%% D + inv(D) = " << endl;
                outfile << mat4 - inverse4;
                outfile << endl;
                outfile << "%% D - D^H = " << endl;
                outfile << mat4 -  transpose4;
                outfile << endl;
                outfile << "%% D * D^H = " << endl;
                outfile << mat4 *  transpose4;
                outfile << endl;
      /*    outfile << "D / inv(D) = " << endl;
                outfile << mat4 /  inverse4;
                outfile << endl;*/
                outfile << "%% D(1,2) = " << endl;
                outfile << mat4(0,1) << endl;
                outfile << endl;
      
                outfile << "%% E = " << endl;
                outfile << mat5;
                outfile << endl;
                outfile << "%% det(E) = " << deter5 << endl;
                outfile << endl;
                outfile << "%% inv(E) = " << endl;
                outfile << inverse5;
                outfile << endl;
                outfile << "%% E^H = " << endl;
                outfile << transpose5;
                outfile << endl;
                outfile << "%% 3*E = " << endl;
                outfile << f2*mat5;
                outfile << endl;
                outfile << "%% E + inv(E) = " << endl;
                outfile << mat5 - inverse5;
                outfile << endl;
                outfile << "%% E - E^H = " << endl;
                outfile << mat5 -  transpose5;
                outfile << endl;
                outfile << "%% E * E^H = " << endl;
                outfile << mat5 *  transpose5;
                outfile << endl;
      /*    outfile << "E / inv(E) = " << endl;
                outfile << mat5 /  inverse5;
                outfile << endl;*/
                outfile << "%% E(1,2) = " << endl;
                outfile << mat5(0,1) << endl;
                outfile << endl;
      
                outfile << "%% F = " << endl;
                outfile << mat6;
                outfile << endl;
                outfile << "%% det(F) = " << deter6 << endl;
                outfile << endl;
                outfile << "%% inv(F) = " << endl;
                outfile << inverse6;
                outfile << endl;
                outfile << "%% F^H = " << endl;
                outfile << transpose6;
                outfile << endl;
                outfile << "%% 3*F = " << endl;
                outfile << f2*mat6;
                outfile << endl;
                outfile << "%% F + inv(F) = " << endl;
                outfile << mat6 - inverse6;
                outfile << endl;
                outfile << "%% F - F^H = " << endl;
                outfile << mat6 -  transpose6;
                outfile << endl;
                outfile << "%% F * F^H = " << endl;
                outfile << mat6 *  transpose6;
                outfile << endl;
      /*    outfile << "F / inv(F) = " << endl;
                outfile << mat6 /  inverse6;
                outfile << endl;*/
                outfile << "%% F(1,2) = " << endl;
                outfile << mat6(0,1) << endl;
                outfile << endl;
      
                outfile << "%% G = " << endl;
                outfile << mat7;
                outfile << endl;
                outfile << "%% det(G) = " <<  real(deter7) << "+" << imag(deter7) << "i" << endl;
                outfile << endl;
                outfile << "%% inv(G) = " << endl;
                outfile << inverse7;
                outfile << endl;
                outfile << "%% G^H = " << endl;
                outfile << transpose7;
                outfile << endl;
                outfile << "%% 3.5*G = " << endl;
                outfile << f3*mat7;
                outfile << endl;
                outfile << "%% G + inv(G) = " << endl;
                outfile << mat7 - inverse7;
                outfile << endl;
                outfile << "%% G - G^H = " << endl;
                outfile << mat7 -  transpose7;
                outfile << endl;
                outfile << "%% G * G^H = " << endl;
                outfile << mat7 *  transpose7;
                outfile << endl;
                outfile << "%% G / inv(G) = " << endl;
                outfile << mat7 /  inverse7;
                outfile << endl;
                outfile << "%% G(1,2) = " << endl;
                outfile << real(mat7(0,1)) << "+" << imag(mat7(0,1)) << "i" << endl;
                outfile << endl;
      
                outfile << "%% H = " << endl;
                outfile << mat8;
                outfile << endl;
                outfile << "%% det(H) = " <<  real(deter8) << "+" << imag(deter8) << "i" << endl;
                outfile << endl;
                outfile << "%% inv(H) = " << endl;
                outfile << inverse8;
                outfile << endl;
                outfile << "%% H^H = " << endl;
                outfile << transpose8;
                outfile << endl;
                outfile << "%% 3.5*H = " << endl;
                outfile << f3*mat8;
                outfile << endl;
                outfile << "%% H + inv(H) = " << endl;
                outfile << mat8 - inverse8;
                outfile << endl;
                outfile << "%% H - H^H = " << endl;
                outfile << mat8 -  transpose8;
                outfile << endl;
                outfile << "%% H * H^H = " << endl;
                outfile << mat8 *  transpose8;
                outfile << endl;
                outfile << "%% H / inv(H) = " << endl;
                outfile << mat8 /  inverse8;
                outfile << endl;
                outfile << "%% H(1,2) = " << endl;
                outfile << real(mat8(0,1)) << "+" << imag(mat8(0,1)) << "i" << endl;
                outfile << endl;
      
                outfile << "%% I = " << endl;
                outfile << mat9;
                outfile << endl;
                outfile << "%% det(I) = " <<  real(deter9) << "+" << imag(deter9) << "i" << endl;
                outfile << endl;
                outfile << "%% inv(I) = " << endl;
                outfile << inverse9;
                outfile << endl;
                outfile << "%% I^H = " << endl;
                outfile << transpose9;
                outfile << endl;
                outfile << "%% 3.5*I = " << endl;
                outfile << f3*mat9;
                outfile << endl;
                outfile << "%% I + inv(I) = " << endl;
                outfile << mat9 - inverse9;
                outfile << endl;
                outfile << "%% I - I^H = " << endl;
                outfile << mat9 -  transpose9;
                outfile << endl;
                outfile << "%% I * I^H = " << endl;
                outfile << mat9 *  transpose9;
                outfile << endl;
                outfile << "%% I / inv(I) = " << endl;
                outfile << mat9 /  inverse9;
                outfile << endl;
                outfile << "%% I(1,2) = " << endl;
                outfile << real(mat9(0,1)) << "+" << imag(mat9(0,1)) << "i" << endl;
                outfile << endl;
                outfile.close();
            }
        }
      
        if (print_csv){
            string OutFileName = "MatrixTest1.csv";
            string OutFilePathString = WorkDir + "/" + OutDir + "/" + OutFileName;
            const char * OutFilePath = OutFilePathString.c_str();
            ofstream outfile;
            outfile.open(OutFilePath);
            if(outfile.is_open()){
                
                //outfile << mat1;
                outfile << deter1 << endl;
                outfile << inverse1;
                outfile << transpose1;
                outfile << f1*mat1;
                outfile << mat1 - inverse1;
                outfile << mat1 -  transpose1;
                outfile << mat1 *  transpose1;
                outfile << mat1 /  inverse1;
                outfile << mat1(0,1) << endl;
      
                //outfile << mat2;
                outfile << deter2 << endl;
                outfile << inverse2;
                outfile << transpose2;
                outfile << f1*mat2;
                outfile << mat2 - inverse2;
                outfile << mat2 -  transpose2;
                outfile << mat2 *  transpose2;
                outfile << mat2 /  inverse2;
                outfile << mat2(0,1) << endl;
                
                //outfile << mat3;
                outfile << deter3 << endl;
                outfile << inverse3;
                outfile << transpose3;
                outfile << f1*mat3;
                outfile << mat3 - inverse3;
                outfile << mat3 -  transpose3;
                outfile << mat3 *  transpose3;
                outfile << mat3 /  inverse3;
                outfile << mat3(0,1) << endl;
                
                //outfile << mat4;
                outfile << deter4 << endl;
                outfile << inverse4;
                outfile << transpose4;
                outfile << f2*mat4;
                outfile << mat4 - inverse4;
                outfile << mat4 -  transpose4;
                outfile << mat4 *  transpose4;
                outfile << mat4;
                outfile << mat4(0,1) << endl;
                
                //outfile << mat5;
                outfile << deter5 << endl;
                outfile << inverse5;
                outfile << transpose5;
                outfile << f2*mat5;
                outfile << mat5 - inverse5;
                outfile << mat5 -  transpose5;
                outfile << mat5 *  transpose5;
                outfile << mat5;
                outfile << mat5(0,1) << endl;
                
                //outfile << mat6;
                outfile << deter6 << endl;
                outfile << inverse6;
                outfile << transpose6;
                outfile << f2*mat6;
                outfile << mat6 - inverse6;
                outfile << mat6 -  transpose6;
                outfile << mat6 *  transpose6;
                outfile << mat6;
                outfile << mat6(0,1) << endl;
                
                //outfile << mat7;
                outfile <<  real(deter7) << "+" << imag(deter7) << "i" << endl;
                outfile << inverse7;
                outfile << transpose7;
                outfile << f3*mat7;
                outfile << mat7 - inverse7;
                outfile << mat7 -  transpose7;
                outfile << mat7 *  transpose7;
                outfile << mat7 /  inverse7;
                outfile << real(mat7(0,1)) << "+" << imag(mat7(0,1)) << "i" << endl;
      
                //outfile << mat8;
                outfile <<  real(deter8) << "+" << imag(deter8) << "i" << endl;
                outfile << inverse8;
                outfile << transpose8;
                outfile << f3*mat8;
                outfile << mat8 - inverse8;
                outfile << mat8 -  transpose8;
                outfile << mat8 *  transpose8;
                outfile << mat8 /  inverse8;
                outfile << real(mat8(0,1)) << "+" << imag(mat8(0,1)) << "i" << endl;
                
                //outfile << mat9;
                outfile <<  real(deter9) << "+" << imag(deter9) << "i" << endl;
                outfile << inverse9;
                outfile << transpose9;
                outfile << f3*mat9;
                outfile << mat9 - inverse9;
                outfile << mat9 -  transpose9;
                outfile << mat9 *  transpose9;
                outfile << mat9 /  inverse9;
                outfile << real(mat9(0,1)) << "+" << imag(mat9(0,1)) << "i" << endl;
                
                outfile.close();
            }
        }
        
        if(print_mat){
            string MFileName = "MatrixTest1.m";
            string InpName = "MatrixTest1.csv";
            string MFilePathString = WorkDir + "/" + OutDir + "/" + MFileName;
            const char * MFilePath = MFilePathString.c_str();
            ofstream mfile;
            mfile.open(MFilePath);
            if(mfile.is_open()){
                mfile << "%%" << endl;
                mfile << "%" << endl;
                mfile << "%********************************************************************************* " << endl;
                mfile << "%                  " << date << " (GMT)" << endl;
                mfile << "%" << endl;
                mfile << "%" << endl;
                mfile << "%                                    EUSMAT" << endl;
                mfile << "%                        European School of Materials" << endl;
                mfile << "%" << endl;
                mfile << "%                                    DocMASE" << endl;
                mfile << "%                 DOCTORATE IN MATERIALS SCIENCE AND ENGINEERING" << endl;
                mfile << "%" << endl;
                mfile << "%" << endl;
                mfile << "%      MECHANICS OF EXTREME THIN COMPOSITE LAYERS FOR AEROSPACE APPLICATIONS" << endl;
                mfile << "%" << endl;
                mfile << "%" << endl;
                mfile << "%                        Differential geometry tools" << endl;
                mfile << "%" << endl;
                mfile << "%                                     by" << endl;
                mfile << "%" << endl;
                mfile << "%                            Luca Di Stasio, 2016" << endl;
                mfile << "%" << endl;
                mfile << "%" << endl;
                mfile << "%********************************************************************************* " << endl;
                mfile << "%" << endl;
                mfile << "% Unit tests for the C++ class matrix.h" << endl;
                mfile << "%" << endl;
                mfile << "%%" << endl;
                mfile << "" << endl;
                mfile << "clear all" << endl;
                mfile << "close all" << endl;
                mfile << "clc" << endl;
                mfile << "" << endl;
                mfile << "" << endl;
                mfile << "inputfilename = '" << InpName << "';" << endl;
                mfile << "" << endl;
                mfile << "f1 = 3.5;" << endl;
                mfile << "f2 = 3;" << endl;
                mfile << "f3 = 3.5+0i;" << endl;
                mfile << "" << endl;
                mfile << "tol = 10^-14;" << endl;
                mfile << "reltol = 10^-15;" << endl;
                mfile << "" << endl;
                mfile << "tol1 = 10^-13;" << endl;
                mfile << "reltol1 = 10^-14;" << endl;
                mfile << "" << endl;
                mfile << "tol2 = 10^-12;" << endl;
                mfile << "reltol2 = 10^-13;" << endl;
                mfile << "" << endl;
                mfile << "tol3 = 10^-11;" << endl;
                mfile << "reltol3 = 10^-12;" << endl;
                mfile << "" << endl;
                mfile << "tol4 = 10^-9;" << endl;
                mfile << "reltol4 = 10^-10;" << endl;
                mfile << "" << endl;
                mfile << "tol5 = 10^-8;" << endl;
                mfile << "reltol5 = 10^-9;" << endl;
                mfile << "" << endl;
                mfile << "tol6 = 10^-7;" << endl;
                mfile << "reltol6 = 10^-8;" << endl;
                mfile << "" << endl;
                mfile << "tol7 = 10^-6;" << endl;
                mfile << "reltol7 = 10^-7;" << endl;
                mfile << "" << endl;
                mfile << "tol8 = 10^-5;" << endl;
                mfile << "reltol8 = 10^-6;" << endl;
                mfile << "" << endl;
                mfile << "tol9 = 10^-4;" << endl;
                mfile << "reltol9 = 10^-5;" << endl;
                mfile << "" << endl;
                mfile << "tol10 = 10^-3;" << endl;
                mfile << "reltol10 = 10^-4;" << endl;
                mfile << "" << endl;
                vector<int> dim;
                for(unsigned int k=0; k<3; k++){
                    dim = doublemat[k].size();
                    mfile << names[k] << "_matlab = [";
                    for(unsigned int i=0; i<dim[0]; i++){
                         for(unsigned int j=0; j<dim[1]; j++){
                            mfile << doublemat[k](i,j);
                            if(j==dim[1]-1 && i==dim[0]-1){
                                mfile << "];" << endl;
                            }else if(j==dim[1]-1 && i!=dim[0]-1){
                                mfile << ";";
                            }else{
                                mfile << ",";
                            }
                        }
                    }
                }
                for(unsigned int k=0; k<3; k++){
                    dim = intmat[k].size();
                    mfile << names[k+3] << "_matlab = [";
                    for(unsigned int i=0; i<dim[0]; i++){
                        for(unsigned int j=0; j<dim[1]; j++){
                            mfile << intmat[k](i,j);
                            if(j==dim[1]-1 && i==dim[0]-1){
                                mfile << "];" << endl;
                            }else if(j==dim[1]-1 && i!=dim[0]-1){
                                mfile << ";";
                            }else{
                                mfile << ",";
                            }
                        }
                    }
                }
                for(unsigned int k=0; k<3; k++){
                    dim = complexmat[k].size();
                    mfile << names[k+6] << "_matlab = [";
                    for(unsigned int i=0; i<dim[0]; i++){
                        for(unsigned int j=0; j<dim[1]; j++){
                            if(imag(complexmat[k](i,j))<0.0){
                                mfile << real(complexmat[k](i,j)) << imag(complexmat[k](i,j)) << "i" ;
                            }else{
                                mfile << real(complexmat[k](i,j)) << "+" << imag(complexmat[k](i,j)) << "i" ;
                            }
                            if(j==dim[1]-1 && i==dim[0]-1){
                                mfile << "];" << endl;
                            }else if(j==dim[1]-1 && i!=dim[0]-1){
                                mfile << ";";
                            }else{
                                mfile << ",";
                            }
                        }
                    }
                }
                mfile << endl;
                for(unsigned int k=0; k<9; k++){
                    mfile << "det_" << names[k] << "_matlab = det(" << names[k] << "_matlab);" << endl;
                    mfile << "inv_" << names[k] << "_matlab = inv(" << names[k] << "_matlab);" << endl;
                    mfile << "transpose_" << names[k] << "_matlab = " << names[k] << "_matlab';" << endl;
                    mfile << "scalarmult_" << names[k] << "_matlab = f" << k/3 + 1 << "*" << names[k] << "_matlab;" << endl;
                    mfile << "op1_" << names[k] << "_matlab = " << names[k] << "_matlab - inv_" << names[k] << "_matlab;" << endl;
                    mfile << "op2_" << names[k] << "_matlab = " << names[k] << "_matlab - transpose_" << names[k] << "_matlab;" << endl;
                    mfile << "op3_" << names[k] << "_matlab = " << names[k] << "_matlab .* transpose_" << names[k] << "_matlab;" << endl;
                    if(k!=3 && k!=4 && k!=5){
                        mfile << "op4_" << names[k] << "_matlab = " << names[k] << "_matlab ./ inv_" << names[k] << "_matlab;" << endl;
                    }
                    mfile << "elementof_" << names[k] << "_matlab = " << names[k] << "_matlab(1,2);" << endl;
                    mfile << endl;
                }
                int R1,C1,R2,C2;
                int sizes[] = {2,3,4,2,3,4,2,3,4};
                int row = 0;
                for(unsigned int k=0; k<9; k++){
                    R1 = row;
                    C1 = 0;
                    R2 = row;
                    C2 = 0;
                    mfile << "det_" << names[k] << "_cpp = csvread(inputfilename," << R1 << "," << C1 << ",[" << R1 << "," << C1 << "," << R2 << "," << C2 << "]);" << endl;
                    R1 = row+0*sizes[k]+1;
                    C1 = 0;
                    R2 = row+1*sizes[k];
                    C2 = 1+k%3;
                    mfile << "inv_" << names[k] << "_cpp = csvread(inputfilename," << R1 << "," << C1 << ",[" << R1 << "," << C1 << "," << R2 << "," << C2 << "]);" << endl;
                    R1 = row+1*sizes[k]+1;
                    C1 = 0;
                    R2 = row+2*sizes[k];
                    C2 = 1+k%3;
                    mfile << "transpose_" << names[k] << "_cpp = csvread(inputfilename," << R1 << "," << C1 << ",[" << R1 << "," << C1 << "," << R2 << "," << C2 << "]);" << endl;
                    R1 = row+2*sizes[k]+1;
                    C1 = 0;
                    R2 = row+3*sizes[k];
                    C2 = 1+k%3;
                    mfile << "scalarmult_" << names[k] << "_cpp = csvread(inputfilename," << R1 << "," << C1 << ",[" << R1 << "," << C1 << "," << R2 << "," << C2 << "]);" << endl;
                    R1 = row+3*sizes[k]+1;
                    C1 = 0;
                    R2 = row+4*sizes[k];
                    C2 = 1+k%3;
                    mfile << "op1_" << names[k] << "_cpp = csvread(inputfilename," << R1 << "," << C1 << ",[" << R1 << "," << C1 << "," << R2 << "," << C2 << "]);" << endl;
                    R1 = row+4*sizes[k]+1;
                    C1 = 0;
                    R2 = row+5*sizes[k];
                    C2 = 1+k%3;
                    mfile << "op2_" << names[k] << "_cpp= csvread(inputfilename," << R1 << "," << C1 << ",[" << R1 << "," << C1 << "," << R2 << "," << C2 << "]);" << endl;
                    R1 = row+5*sizes[k]+1;
                    C1 = 0;
                    R2 = row+6*sizes[k];
                    C2 = 1+k%3;
                    mfile << "op3_" << names[k] << "_cpp = csvread(inputfilename," << R1 << "," << C1 << ",[" << R1 << "," << C1 << "," << R2 << "," << C2 << "]);" << endl;
                    if(k!=3 && k!=4 && k!=5){
                        R1 = row+6*sizes[k]+1;
                        C1 = 0;
                        R2 = row+7*sizes[k];
                        C2 = 1+k%3;
                        mfile << "op4_" << names[k] << "_cpp = csvread(inputfilename," << R1 << "," << C1 << ",[" << R1 << "," << C1 << "," << R2 << "," << C2 << "]);" << endl;
                    }
                    R1 = row+7*sizes[k]+1;
                    C1 = 0;
                    R2 = row+7*sizes[k]+1;
                    C2 = 0;
                    mfile << "elementof_" << names[k] << "_cpp= csvread(inputfilename," << R1 << "," << C1 << ",[" << R1 << "," << C1 << "," << R2 << "," << C2 << "]);" << endl;
                    mfile << endl;
                    row += 2 + 7*sizes[k];
                }
                for(unsigned int k=0; k<9; k++){
                    mfile << "det_" << names[k] << "_diff = " << "det_" << names[k] << "_matlab" << " - " << "det_" << names[k] << "_cpp;" << endl;
                    mfile << "inv_" << names[k] << "_diff = " << "inv_" << names[k] << "_matlab" << " - " << "inv_" << names[k] << "_cpp;" << endl;
                    mfile << "transpose_" << names[k] << "_diff = " << "transpose_" << names[k] << "_matlab" << " - " << "transpose_" << names[k] << "_cpp;" << endl;
                    mfile << "scalarmult_" << names[k] << "_diff = " << "scalarmult_" << names[k] << "_matlab" << " - " << "scalarmult_" << names[k] << "_cpp;" << endl;
                    mfile << "op1_" << names[k] << "_diff = " << "op1_" << names[k] << "_matlab" << " - " << "op1_" << names[k] << "_cpp;" << endl;
                    mfile << "op2_" << names[k] << "_diff = " << "op2_" << names[k] << "_matlab" << " - " << "op2_" << names[k] << "_cpp;" << endl;
                    mfile << "op3_" << names[k] << "_diff = " << "op3_" << names[k] << "_matlab" << " - " << "op3_" << names[k] << "_cpp;" << endl;
                    if(k!=3 && k!=4 && k!=5){
                        mfile << "op4_" << names[k] << "_diff = " << "op4_" << names[k] << "_matlab" << " - " << "op4_" << names[k] << "_cpp;" << endl;
                    }
                    mfile << "elementof_" << names[k] << "_diff = " << "elementof_" << names[k] << "_matlab" << " - " << "elementof_" << names[k] << "_cpp;" << endl;
                    mfile << endl;
                }
                for(unsigned int k=0; k<9; k++){
                    mfile << "det_" << names[k] << "_reldiff = " << "det_" << names[k] << "_diff" << " ./ " << "det_" << names[k] << "_matlab;" << endl;
                    mfile << "inv_" << names[k] << "_reldiff = " << "inv_" << names[k] << "_diff" << " ./ " << "inv_" << names[k] << "_matlab;" << endl;
                    mfile << "transpose_" << names[k] << "_reldiff = " << "transpose_" << names[k] << "_diff" << " ./ " << "transpose_" << names[k] << "_matlab;" << endl;
                    mfile << "scalarmult_" << names[k] << "_reldiff = " << "scalarmult_" << names[k] << "_diff" << " ./ " << "scalarmult_" << names[k] << "_matlab;" << endl;
                    mfile << "op1_" << names[k] << "_reldiff = " << "op1_" << names[k] << "_diff" << " ./ " << "op1_" << names[k] << "_matlab;" << endl;
                    mfile << "op2_" << names[k] << "_reldiff = " << "op2_" << names[k] << "_diff" << " ./ " << "op2_" << names[k] << "_matlab;" << endl;
                    mfile << "op3_" << names[k] << "_reldiff = " << "op3_" << names[k] << "_diff" << " ./ " << "op3_" << names[k] << "_matlab;" << endl;
                    if(k!=3 && k!=4 && k!=5){
                        mfile << "op4_" << names[k] << "_reldiff = " << "op4_" << names[k] << "_diff" << " ./ " << "op4_" << names[k] << "_matlab;" << endl;
                    }
                    mfile << "elementof_" << names[k] << "_reldiff = " << "elementof_" << names[k] << "_diff" << " ./ " << "elementof_" << names[k] << "_matlab;" << endl;
                    mfile << endl;
                }
                for(unsigned int k=0; k<9; k++){
                    mfile << "det_" << names[k] << "_normdiff = " << "norm(det_" << names[k] << "_diff);" << endl;
                    mfile << "inv_" << names[k] << "_normdiff = " << "norm(inv_" << names[k] << "_diff);" << endl;
                    mfile << "transpose_" << names[k] << "_normdiff = " << "norm(transpose_" << names[k] << "_diff);" << endl;
                    mfile << "scalarmult_" << names[k] << "_normdiff = " << "norm(scalarmult_" << names[k] << "_diff);" << endl;
                    mfile << "op1_" << names[k] << "_normdiff = " << "norm(op1_" << names[k] << "_diff);" << endl;
                    mfile << "op2_" << names[k] << "_normdiff = " << "norm(op2_" << names[k] << "_diff);" << endl;
                    mfile << "op3_" << names[k] << "_normdiff = " << "norm(op3_" << names[k] << "_diff);" << endl;
                    if(k!=3 && k!=4 && k!=5){
                        mfile << "op4_" << names[k] << "_normdiff = " << "norm(op4_" << names[k] << "_diff);" << endl;
                    }
                    mfile << "elementof_" << names[k] << "_normdiff = " << "norm(elementof_" << names[k] << "_diff);" << endl;
                    mfile << endl;
                }
                for(unsigned int k=0; k<9; k++){
                    mfile << "det_" << names[k] << "_normreldiff = " << "norm(det_" << names[k] << "_reldiff);" << endl;
                    mfile << "inv_" << names[k] << "_normreldiff = " << "norm(inv_" << names[k] << "_reldiff);" << endl;
                    mfile << "transpose_" << names[k] << "_normreldiff = " << "norm(transpose_" << names[k] << "_reldiff);" << endl;
                    mfile << "scalarmult_" << names[k] << "_normreldiff = " << "norm(scalarmult_" << names[k] << "_reldiff);" << endl;
                    mfile << "op1_" << names[k] << "_normreldiff = " << "norm(op1_" << names[k] << "_reldiff);" << endl;
                    mfile << "op2_" << names[k] << "_normreldiff = " << "norm(op2_" << names[k] << "_reldiff);" << endl;
                    mfile << "op3_" << names[k] << "_normreldiff = " << "norm(op3_" << names[k] << "_reldiff);" << endl;
                    if(k!=3 && k!=4 && k!=5){
                        mfile << "op4_" << names[k] << "_normreldiff = " << "norm(op4_" << names[k] << "_reldiff);" << endl;
                    }
                    mfile << "elementof_" << names[k] << "_normreldiff = " << "norm(elementof_" << names[k] << "_reldiff);" << endl;
                    mfile << endl;
                }
                for(unsigned int k=0; k<9; k++){
                    mfile << "" << endl;
                    mfile << "disp('')" << endl;
                    mfile << "disp('=========================== MATRIX "<< k+1 << "===========================')" << endl;
                    mfile << "disp('')" << endl;
                    mfile << "" << endl;
                    mfile << "if(det_" << names[k] << "_normreldiff" << "<=reltol && det_" << names[k] << "_normdiff" << "<=tol)" << endl;
                    mfile << "  disp('Test " << k*9+1 << " (determinant): OK  ----> rel_tol < 10^-15   &   tol < 10^-14')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(det_" << names[k] << "_normreldiff" << "<=reltol" << m << " && det_" << names[k] << "_normdiff" << "<=tol" << m << ")" << endl;
                        mfile << "  disp('Test " << k*9+1 << " (determinant): OK  ----> rel_tol < 10^-" << 15-m << "   &   tol < 10^-" << 14-m << "')" << endl;
                    }
                    mfile << "elseif(det_" << names[k] << "_normreldiff" << "<=reltol)" << endl;
                    mfile << "  disp('Test " << k*9+1 << " (determinant): OK  ----> rel_tol < 10^-15')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(det_" << names[k] << "_normreldiff" << "<=reltol" << m << ")" << endl;
                        mfile << "  disp('Test " << k*9+1 << " (determinant): OK  ----> rel_tol < 10^-" << 15-m << "')" << endl;
                    }
                    mfile << "elseif(det_" << names[k] << "_normdiff" << "<=tol && isnan(det_" << names[k] << "_normreldiff))" << endl;
                    mfile << "  disp('Test " << k*9+1 << " (determinant): OK  ----> tol < 10^-15    (NOTE: normreldiff is NaN, probably due to division by zero)')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(det_" << names[k] << "_normdiff" << "<=tol" << m << " && isnan(det_" << names[k] << "_normreldiff))" << endl;
                        mfile << "  disp('Test " << k*9+1 << " (determinant): OK  ----> rel_tol < 10^-" << 15-m << "    (NOTE: normreldiff is NaN, probably due to division by zero)')" << endl;
                    }
                    mfile << "else" << endl;
                    mfile << "  disp('Test " << k*9+1 << " (determinant): =================> !!WARNING!! <=================')" << endl;
                    mfile << "end" << endl;
                    mfile << "" << endl;
    //================================================================================================================================================================================================//
                    mfile << "if(inv_" << names[k] << "_normreldiff" << "<=reltol && inv_" << names[k] << "_normdiff" << "<=tol)" << endl;
                    mfile << "  disp('Test " << k*9+2 << " (inverse): OK  ----> rel_tol < 10^-15   &   tol < 10^-14')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(inv_" << names[k] << "_normreldiff" << "<=reltol" << m << " && inv_" << names[k] << "_normdiff" << "<=tol" << m << ")" << endl;
                        mfile << "  disp('Test " << k*9+2 << " (inverse): OK  ----> rel_tol < 10^-" << 15-m << "   &   tol < 10^-" << 14-m << "')" << endl;
                    }
                    mfile << "elseif(inv_" << names[k] << "_normreldiff" << "<=reltol)" << endl;
                    mfile << "  disp('Test " << k*9+2 << " (inverse): OK  ----> rel_tol < 10^-15')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(inv_" << names[k] << "_normreldiff" << "<=reltol" << m << ")" << endl;
                        mfile << "  disp('Test " << k*9+2 << " (inverse): OK  ----> rel_tol < 10^-" << 15-m << "')" << endl;
                    }
                    mfile << "elseif(inv_" << names[k] << "_normdiff" << "<=tol && isnan(inv_" << names[k] << "_normreldiff))" << endl;
                    mfile << "  disp('Test " << k*9+2 << " (inverse): OK  ----> tol < 10^-15    (NOTE: normreldiff is NaN, probably due to division by zero)')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(inv_" << names[k] << "_normdiff" << "<=tol" << m << " && isnan(inv_" << names[k] << "_normreldiff))" << endl;
                        mfile << "  disp('Test " << k*9+2 << " (inverse): OK  ----> rel_tol < 10^-" << 15-m << "    (NOTE: normreldiff is NaN, probably due to division by zero)')" << endl;
                    }
                    mfile << "else" << endl;
                    mfile << "  disp('Test " << k*9+2 << " (inverse): =================> !!WARNING!! <=================')" << endl;
                    mfile << "end" << endl;
                    mfile << "" << endl;
    //================================================================================================================================================================================================//                
                    mfile << "if(transpose_" << names[k] << "_normreldiff" << "<=reltol && transpose_" << names[k] << "_normdiff" << "<=tol)" << endl;
                    mfile << "  disp('Test " << k*9+3 << " (hermitian): OK  ----> rel_tol < 10^-15   &   tol < 10^-14')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(transpose_" << names[k] << "_normreldiff" << "<=reltol" << m << " && transpose_" << names[k] << "_normdiff" << "<=tol" << m << ")" << endl;
                        mfile << "  disp('Test " << k*9+3 << " (hermitian): OK  ----> rel_tol < 10^-" << 15-m << "   &   tol < 10^-" << 14-m << "')" << endl;
                    }
                    mfile << "elseif(transpose_" << names[k] << "_normreldiff" << "<=reltol)" << endl;
                    mfile << "  disp('Test " << k*9+3 << " (hermitian): OK  ----> rel_tol < 10^-15')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(transpose_" << names[k] << "_normreldiff" << "<=reltol" << m << ")" << endl;
                        mfile << "  disp('Test " << k*9+3 << " (hermitian): OK  ----> rel_tol < 10^-" << 15-m << "')" << endl;
                    }
                    mfile << "elseif(transpose_" << names[k] << "_normdiff" << "<=tol && isnan(transpose_" << names[k] << "_normreldiff))" << endl;
                    mfile << "  disp('Test " << k*9+3 << " (hermitian): OK  ----> tol < 10^-15    (NOTE: normreldiff is NaN, probably due to division by zero)')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(transpose_" << names[k] << "_normdiff" << "<=tol" << m << " && isnan(transpose_" << names[k] << "_normreldiff))" << endl;
                        mfile << "  disp('Test " << k*9+3 << " (hermitian): OK  ----> rel_tol < 10^-" << 15-m << "    (NOTE: normreldiff is NaN, probably due to division by zero)')" << endl;
                    }
                    mfile << "else" << endl;
                    mfile << "  disp('Test " << k*9+3 << " (hermitian): =================> !!WARNING!! <=================')" << endl;
                    mfile << "end" << endl;
                    mfile << "" << endl;
    //================================================================================================================================================================================================//                
                    mfile << "if(scalarmult_" << names[k] << "_normreldiff" << "<=reltol && scalarmult_" << names[k] << "_normdiff" << "<=tol)" << endl;
                    mfile << "  disp('Test " << k*9+4 << " (scalar multiplication): OK  ----> rel_tol < 10^-15   &   tol < 10^-14')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(scalarmult_" << names[k] << "_normreldiff" << "<=reltol" << m << " && scalarmult_" << names[k] << "_normdiff" << "<=tol" << m << ")" << endl;
                        mfile << "  disp('Test " << k*9+4 << " (scalar multiplication): OK  ----> rel_tol < 10^-" << 15-m << "   &   tol < 10^-" << 14-m << "')" << endl;
                    }
                    mfile << "elseif(scalarmult_" << names[k] << "_normreldiff" << "<=reltol)" << endl;
                    mfile << "  disp('Test " << k*9+4 << " (scalar multiplication): OK  ----> rel_tol < 10^-15')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(scalarmult_" << names[k] << "_normreldiff" << "<=reltol" << m << ")" << endl;
                        mfile << "  disp('Test " << k*9+4 << " (scalar multiplication): OK  ----> rel_tol < 10^-" << 15-m << "')" << endl;
                    }
                    mfile << "elseif(scalarmult_" << names[k] << "_normdiff" << "<=tol && isnan(scalarmult_" << names[k] << "_normreldiff))" << endl;
                    mfile << "  disp('Test " << k*9+4 << " (scalar multiplication): OK  ----> tol < 10^-15    (NOTE: normreldiff is NaN, probably due to division by zero)')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(scalarmult_" << names[k] << "_normdiff" << "<=tol" << m << " && isnan(scalarmult_" << names[k] << "_normreldiff))" << endl;
                        mfile << "  disp('Test " << k*9+4 << " (scalar multiplication): OK  ----> rel_tol < 10^-" << 15-m << "    (NOTE: normreldiff is NaN, probably due to division by zero)')" << endl;
                    }
                    mfile << "else" << endl;
                    mfile << "  disp('Test " << k*9+4 << " (scalar multiplication): =================> !!WARNING!! <=================')" << endl;
                    mfile << "end" << endl;
                    mfile << "" << endl;
    //================================================================================================================================================================================================//                
                    mfile << "if(op1_" << names[k] << "_normreldiff" << "<=reltol && op1_" << names[k] << "_normdiff" << "<=tol)" << endl;
                    mfile << "  disp('Test " << k*9+5 << " (A + inv(A)): OK  ----> rel_tol < 10^-15   &   tol < 10^-14')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(op1_" << names[k] << "_normreldiff" << "<=reltol" << m << " && op1_" << names[k] << "_normdiff" << "<=tol" << m << ")" << endl;
                        mfile << "  disp('Test " << k*9+5 << " (A + inv(A)): OK  ----> rel_tol < 10^-" << 15-m << "   &   tol < 10^-" << 14-m << "')" << endl;
                    }
                    mfile << "elseif(op1_" << names[k] << "_normreldiff" << "<=reltol)" << endl;
                    mfile << "  disp('Test " << k*9+5 << " (A + inv(A)): OK  ----> rel_tol < 10^-15')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(op1_" << names[k] << "_normreldiff" << "<=reltol" << m << ")" << endl;
                        mfile << "  disp('Test " << k*9+5 << " (A + inv(A)): OK  ----> rel_tol < 10^-" << 15-m << "')" << endl;
                    }
                    mfile << "elseif(op1_" << names[k] << "_normdiff" << "<=tol && isnan(op1_" << names[k] << "_normreldiff))" << endl;
                    mfile << "  disp('Test " << k*9+5 << " (A + inv(A)): OK  ----> tol < 10^-15    (NOTE: normreldiff is NaN, probably due to division by zero)')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(op1_" << names[k] << "_normdiff" << "<=tol" << m << " && isnan(op1_" << names[k] << "_normreldiff))" << endl;
                        mfile << "  disp('Test " << k*9+5 << " (A + inv(A)): OK  ----> rel_tol < 10^-" << 15-m << "    (NOTE: normreldiff is NaN, probably due to division by zero)')" << endl;
                    }
                    mfile << "else" << endl;
                    mfile << "  disp('Test " << k*9+5 << " (A + inv(A)): =================> !!WARNING!! <=================')" << endl;
                    mfile << "end" << endl;
                    mfile << "" << endl;
    //================================================================================================================================================================================================//                
                    mfile << "if(op2_" << names[k] << "_normreldiff" << "<=reltol && op2_" << names[k] << "_normdiff" << "<=tol)" << endl;
                    mfile << "  disp('Test " << k*9+6 << " (A - A^H): OK  ----> rel_tol < 10^-15   &   tol < 10^-14')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(op2_" << names[k] << "_normreldiff" << "<=reltol" << m << " && op2_" << names[k] << "_normdiff" << "<=tol" << m << ")" << endl;
                        mfile << "  disp('Test " << k*9+6 << " (A - A^H): OK  ----> rel_tol < 10^-" << 15-m << "   &   tol < 10^-" << 14-m << "')" << endl;
                    }
                    mfile << "elseif(op2_" << names[k] << "_normreldiff" << "<=reltol)" << endl;
                    mfile << "  disp('Test " << k*9+6 << " (A - A^H): OK  ----> rel_tol < 10^-15')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(op2_" << names[k] << "_normreldiff" << "<=reltol" << m << ")" << endl;
                        mfile << "  disp('Test " << k*9+6 << " (A - A^H): OK  ----> rel_tol < 10^-" << 15-m << "')" << endl;
                    }
                    mfile << "elseif(op2_" << names[k] << "_normdiff" << "<=tol && isnan(op2_" << names[k] << "_normreldiff))" << endl;
                    mfile << "  disp('Test " << k*9+6 << " (A - A^H): OK  ----> tol < 10^-15    (NOTE: normreldiff is NaN, probably due to division by zero)')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(op2_" << names[k] << "_normdiff" << "<=tol" << m << " && isnan(op2_" << names[k] << "_normreldiff))" << endl;
                        mfile << "  disp('Test " << k*9+6 << " (A - A^H): OK  ----> rel_tol < 10^-" << 15-m << "    (NOTE: normreldiff is NaN, probably due to division by zero)')" << endl;
                    }
                    mfile << "else" << endl;
                    mfile << "  disp('Test " << k*9+6 << " (A - A^H): =================> !!WARNING!! <=================')" << endl;
                    mfile << "end" << endl;
                    mfile << "" << endl;
    //================================================================================================================================================================================================//                 
                    mfile << "if(op3_" << names[k] << "_normreldiff" << "<=reltol && op3_" << names[k] << "_normdiff" << "<=tol)" << endl;
                    mfile << "  disp('Test " << k*9+7 << " (A * A^H): OK  ----> rel_tol < 10^-15   &   tol < 10^-14')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(op3_" << names[k] << "_normreldiff" << "<=reltol" << m << " && op3_" << names[k] << "_normdiff" << "<=tol" << m << ")" << endl;
                        mfile << "  disp('Test " << k*9+7 << " (A * A^H): OK  ----> rel_tol < 10^-" << 15-m << "   &   tol < 10^-" << 14-m << "')" << endl;
                    }
                    mfile << "elseif(op3_" << names[k] << "_normreldiff" << "<=reltol)" << endl;
                    mfile << "  disp('Test " << k*9+7 << " (A * A^H): OK  ----> rel_tol < 10^-15')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(op3_" << names[k] << "_normreldiff" << "<=reltol" << m << ")" << endl;
                        mfile << "  disp('Test " << k*9+7 << " (A * A^H): OK  ----> rel_tol < 10^-" << 15-m << "')" << endl;
                    }
                    mfile << "elseif(op3_" << names[k] << "_normdiff" << "<=tol && isnan(op3_" << names[k] << "_normreldiff))" << endl;
                    mfile << "  disp('Test " << k*9+7 << " (A * A^H): OK  ----> tol < 10^-15    (NOTE: normreldiff is NaN, probably due to division by zero)')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(op3_" << names[k] << "_normdiff" << "<=tol" << m << " && isnan(op3_" << names[k] << "_normreldiff))" << endl;
                        mfile << "  disp('Test " << k*9+7 << " (A * A^H): OK  ----> rel_tol < 10^-" << 15-m << "    (NOTE: normreldiff is NaN, probably due to division by zero)')" << endl;
                    }
                    mfile << "else" << endl;
                    mfile << "  disp('Test " << k*9+7 << " (A * A^H): =================> !!WARNING!! <=================')" << endl;
                    mfile << "end" << endl;
                    mfile << "" << endl;
    //================================================================================================================================================================================================//                
                    if(k!=3 && k!=4 && k!=5){
                        mfile << "if(op4_" << names[k] << "_normreldiff" << "<=reltol && op4_" << names[k] << "_normdiff" << "<=tol)" << endl;
                        mfile << "  disp('Test " << k*9+8 << " (A / inv(A)): OK  ----> rel_tol < 10^-15   &   tol < 10^-14')" << endl;
                        for(unsigned int m=1; m<11; m++){
                            mfile << "elseif(op4_" << names[k] << "_normreldiff" << "<=reltol" << m << " && op4_" << names[k] << "_normdiff" << "<=tol" << m << ")" << endl;
                            mfile << "  disp('Test " << k*9+8 << " (A / inv(A)): OK  ----> rel_tol < 10^-" << 15-m << "   &   tol < 10^-" << 14-m << "')" << endl;
                        }
                        mfile << "elseif(op4_" << names[k] << "_normreldiff" << "<=reltol)" << endl;
                        mfile << "  disp('Test " << k*9+8 << " (A / inv(A)): OK  ----> rel_tol < 10^-15')" << endl;
                        for(unsigned int m=1; m<11; m++){
                            mfile << "elseif(op4_" << names[k] << "_normreldiff" << "<=reltol" << m << ")" << endl;
                            mfile << "  disp('Test " << k*9+8 << " (A / inv(A)): OK  ----> rel_tol < 10^-" << 15-m << "')" << endl;
                        }
                        mfile << "elseif(op4_" << names[k] << "_normdiff" << "<=tol && isnan(op4_" << names[k] << "_normreldiff))" << endl;
                        mfile << "  disp('Test " << k*9+8 << " (A / inv(A)): OK  ----> tol < 10^-15    (NOTE: normreldiff is NaN, probably due to division by zero)')" << endl;
                        for(unsigned int m=1; m<11; m++){
                            mfile << "elseif(op4_" << names[k] << "_normdiff" << "<=tol" << m << " && isnan(op4_" << names[k] << "_normreldiff))" << endl;
                            mfile << "  disp('Test " << k*9+8 << " (A / inv(A)): OK  ----> rel_tol < 10^-" << 15-m << "    (NOTE: normreldiff is NaN, probably due to division by zero)')" << endl;
                        }
                        mfile << "else" << endl;
                        mfile << "  disp('Test " << k*9+8 << " (A / inv(A)): =================> !!WARNING!! <=================')" << endl;
                        mfile << "end" << endl;
                    }else{
                        mfile << "  disp('Test not applicable')" << endl;
                    }
                    mfile << "" << endl;
    //================================================================================================================================================================================================//                
                    mfile << "if(elementof_" << names[k] << "_normreldiff" << "<=reltol && elementof_" << names[k] << "_normdiff" << "<=tol)" << endl;
                    mfile << "  disp('Test " << k*9+9 << " (element of): OK  ----> rel_tol < 10^-15   &   tol < 10^-14')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(elementof_" << names[k] << "_normreldiff" << "<=reltol" << m << " && elementof_" << names[k] << "_normdiff" << "<=tol" << m << ")" << endl;
                        mfile << "  disp('Test " << k*9+9 << " (element of): OK  ----> rel_tol < 10^-" << 15-m << "   &   tol < 10^-" << 14-m << "')" << endl;
                    }
                    mfile << "elseif(elementof_" << names[k] << "_normreldiff" << "<=reltol)" << endl;
                    mfile << "  disp('Test " << k*9+9 << " (element of): OK  ----> rel_tol < 10^-15')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(elementof_" << names[k] << "_normreldiff" << "<=reltol" << m << ")" << endl;
                        mfile << "  disp('Test " << k*9+9 << " (element of): OK  ----> rel_tol < 10^-" << 15-m << "')" << endl;
                    }
                    mfile << "elseif(elementof_" << names[k] << "_normdiff" << "<=tol && isnan(elementof_" << names[k] << "_normreldiff))" << endl;
                    mfile << "  disp('Test " << k*9+9 << " (element of): OK  ----> tol < 10^-15    (NOTE: normreldiff is NaN, probably due to division by zero)')" << endl;
                    for(unsigned int m=1; m<11; m++){
                        mfile << "elseif(elementof_" << names[k] << "_normdiff" << "<=tol" << m << " && isnan(elementof_" << names[k] << "_normreldiff))" << endl;
                        mfile << "  disp('Test " << k*9+9 << " (element of): OK  ----> rel_tol < 10^-" << 15-m << "    (NOTE: normreldiff is NaN, probably due to division by zero)')" << endl;
                    }
                    mfile << "else" << endl;
                    mfile << "  disp('Test " << k*9+9 << " (element of): =================> !!WARNING!! <=================')" << endl;
                    mfile << "end" << endl;
    //================================================================================================================================================================================================//                
                    mfile << endl;
                }
                
                mfile.close();
            }
        }
  
    }
    
    if(test2){
        matrix<double> A(3,3);
        matrix<double> B(4,4);
        matrix<double> C(5,5);
        
        matrix<int> D(4,4);
        matrix<int> E(5,5);
        matrix<int> F(6,6);
        
        matrix<complex<double>> G(2,2);
        matrix<complex<double>> H(3,3);
        matrix<complex<double>> J(7,7);
        
        matrix<double> K(3,3);
        matrix<double> I(4,4);
        matrix<double> L(5,5);
        
        matrix<int> M(4,4);
        matrix<int> N(5,5);
        matrix<int> O(6,6);
        
        matrix<complex<double>> P(2,2);
        matrix<complex<double>> Q(3,3);
        matrix<complex<double>> R(7,7);
        
        matrix<double> a_col(3,1);
        matrix<double> b_col(4,1);
        matrix<double> c_col(5,1);
        
        matrix<int> d_col(4,1);
        matrix<int> e_col(5,1);
        matrix<int> f_col(6,1);
        
        matrix<complex<double>> g_col(2,1);
        matrix<complex<double>> h_col(3,1);
        matrix<complex<double>> j_col(7,1);
        
        matrix<double> a_row(1,3);
        matrix<double> b_row(1,4);
        matrix<double> c_row(1,5);
        
        matrix<int> d_row(1,4);
        matrix<int> e_row(1,5);
        matrix<int> f_row(1,6);
        
        matrix<complex<double>> g_row(1,2);
        matrix<complex<double>> h_row(1,3);
        matrix<complex<double>> j_row(1,7);
        
        A.init_random();
        B.init_random(-10,10);
        C.init_random(0.0,100.0);
        D.init_random(-100,100);
        E.init_random(-50.0,50.0);
        F.init_random(-45,5);
        G.init_random();
        H.init_random(-100,100);
        J.init_random(20.5,75.5);
        
        int diag1[] = {5,8,6,-4,-7};
        complex<double> diag2[] = {{1.0,2.5},{7.8,6.1},{7.4,8.7}};
        complex<double> diag3[] = {{-9.0,4.5},{-7.8,18.1},{-7.4,-8.7}};
        
        K.identity();
        I.identity(6.67);
        L.ones();
        M.ones(8.5);
        N.identity();
        O.diag(diag1,5,-1);
        P.identity(complex<double>(8.1,-5.4));
        Q.diag(diag2,3);
        R.diag(diag3,3,4);
        
        a_col.init_random();
        b_col.init_random(-10,10);
        c_col.init_random(0.0,100.0);
        d_col.init_random(-100,100);
        e_col.init_random(-50.0,50.0);
        f_col.init_random(-45,5);
        g_col.init_random();
        h_col.init_random(-100,100);
        j_col.init_random(20.5,75.5);
        
        a_row.init_random();
        b_row.init_random(-10,10);
        c_row.init_random(0.0,100.0);
        d_row.init_random(-100,100);
        e_row.init_random(-50.0,50.0);
        f_row.init_random(-45,5);
        g_row.init_random();
        h_row.init_random(-100,100);
        j_row.init_random(20.5,75.5);
        
        string name;
        
        if(print_screen){
            cout << A.matlab_writer(names[0]);
            cout << endl;
            cout << B.matlab_writer(names[1]);
            cout << endl;
            cout << C.matlab_writer(names[2]);
            cout << endl;
            cout << D.matlab_writer(names[3]);
            cout << endl;
            cout << E.matlab_writer(names[4]);
            cout << endl;
            cout << F.matlab_writer(names[5]);
            cout << endl;
            cout << G.matlab_writer(names[6]);
            cout << endl;
            cout << H.matlab_writer(names[7]);
            cout << endl;
            cout << J.matlab_writer(names[8]);
            cout << endl;
            
            name = names[0] + "_col";
            cout << a_col.matlab_writer(name);
            cout << endl;
            name = names[1] + "_col";
            cout << b_col.matlab_writer(name);
            cout << endl;
            name = names[2] + "_col";
            cout << c_col.matlab_writer(name);
            cout << endl;
            name = names[3] + "_col";
            cout << d_col.matlab_writer(name);
            cout << endl;
            name = names[4] + "_col";
            cout << e_col.matlab_writer(name);
            cout << endl;
            name = names[5] + "_col";
            cout << f_col.matlab_writer(name);
            cout << endl;
            name = names[6] + "_col";
            cout << g_col.matlab_writer(name);
            cout << endl;
            name = names[7] + "_col";
            cout << h_col.matlab_writer(name);
            cout << endl;
            name = names[8] + "_col";
            cout << j_col.matlab_writer(name);
            cout << endl;
            
            name = names[0] + "_row";
            cout << a_row.matlab_writer(name);
            cout << endl;
            name = names[1] + "_row";
            cout << b_row.matlab_writer(name);
            cout << endl;
            name = names[2] + "_row";
            cout << c_row.matlab_writer(name);
            cout << endl;
            name = names[3] + "_row";
            cout << d_row.matlab_writer(name);
            cout << endl;
            name = names[4] + "_row";
            cout << e_row.matlab_writer(name);
            cout << endl;
            name = names[5] + "_row";
            cout << f_row.matlab_writer(name);
            cout << endl;
            name = names[6] + "_row";
            cout << g_row.matlab_writer(name);
            cout << endl;
            name = names[7] + "_row";
            cout << h_row.matlab_writer(name);
            cout << endl;
            name = names[8] + "_row";
            cout << j_row.matlab_writer(name);
            cout << endl;
            
            name = names[9];
            cout << K.matlab_writer(name);
            cout << endl;
            name = names[10];
            cout << I.matlab_writer(name);
            cout << endl;
            name = names[11];
            cout << L.matlab_writer(name);
            cout << endl;
            name = names[12];
            cout << M.matlab_writer(name);
            cout << endl;
            name = names[13];
            cout << N.matlab_writer(name);
            cout << endl;
            name = names[14];
            cout << O.matlab_writer(name);
            cout << endl;
            name = names[15];
            cout << P.matlab_writer(name);
            cout << endl;
            name = names[16];
            cout << Q.matlab_writer(name);
            cout << endl;
            name = names[17];
            cout << R.matlab_writer(name);
            cout << endl;
        }
 
    }
    
    return 0;
}
#ifndef TC_FACTORIZACION_H
#define TC_FACTORIZACION_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <stdio.h>
#include <sstream>

using namespace std;

vector<double> sumaVec(const vector<double> &vec1, const vector<double> &vec2);
vector<double> restaVec(const vector<double> &vec1, const vector<double> &vec2);
double norma1(const vector<double> &v);
void normalizar2(vector<double>& v);
double producto_interno(const vector<double> &v, const vector<double> &v2);
vector<double> mult_matr_por_vect(const vector<vector<double> > &M, const vector<double> &v);
pair<double,vector<double> > metodoPotencia(const vector<vector<double> > &M);
void multVecEsc(vector<double> &vec, double escalar);
void multMatEsc(vector<vector<double> > &mat, double escalar);
vector<vector<double> > multVec(const vector<double> &vec1);
vector<vector<double> > calcularXtX (const vector<vector<double> >& X);
void sumMat(vector<vector<double> > &mat1, const vector<vector<double> > &mat2);
vector< pair<double,vector<double> > > deflacion(vector<vector<double> > mat);
void calcular_svd(const vector<vector<double> > &mat,
                        vector<vector<double> > &Ut,    //
                        vector<double> &Sigm,
                        vector<vector<double> > &Vt,
                        const double tolerance,
                        const bool debug);

#endif //TC_FACTORIZACION_H

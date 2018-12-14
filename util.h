#ifndef TC_UTIL_H
#define TC_UTIL_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <stdio.h>
#include <sstream>
#include "dirent.h"
#include "VectorMapMatrix.h"
#include "calcular_rayos.h"
#include "factorizacion.h"

typedef unsigned char uchar;

#define PI 3.1415926536

const ulong MAX_u_cuadrado = pow(255, 2);


using namespace std;

vector<vector<double>>* leerCSV(string nombreArchivo);
vector<vector<double>> discretizar(const vector<vector<double> >& mat, uint val);
map<uint, double> pasarAMap(const vector<vector<double> >& mat);
vector<double> pasarAVector(const vector<vector<double> >& mat);
VectorMapMatrix  generarRayos(size_t tamMatriz, bool fijos);
VectorMapMatrix  generarRayos_barrido_H(size_t tamMatriz, size_t cada_cuanto);
vector<double> uniformNoise(const vector<double>& t, double init, double end, double sign);
VectorMapMatrix getTraspuesta(const VectorMapMatrix &W);
vector<vector<double> > calcularXtX (const vector<vector<double> >& X);
vector<double> CML(vector<vector<double>> &mat, vector<double> bb,const double tolerancia,const bool debug);
double operator*(const vector<double>& u, const vector<double>& v);
vector<double> operator*(const vector<vector<double> >& M, const vector<double>& v);
void listarDirectorio(const string& directorio,  vector<string>& v);
void escribirVector(string nombreArchivo, vector<double>& vector);
void escribirVectorDeVectores(string nombreArchivo, vector<vector<double>>& vector);
void escribirCSV(string nombreArchivo, vector<double>& vector, size_t ancho);


vector<double> AWGNNoise(const vector<double>& t, const vector<double>& imagen, double porcentajeDeRuido);
vector<double> MWGNNoise(const vector<double>& t, const vector<double>& imagen, double porcentajeDeRuido);

#endif //TC_UTIL_H

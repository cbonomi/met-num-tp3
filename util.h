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
#include "rdtsc.h"
#include "ppmloader.h"
#include "VectorMapMatrix.h"
#include "calcular_rayos.h"

using namespace std;

vector<vector<double>>* leerCSV(string nombreArchivo);
vector<vector<double>> discretizar(const vector<vector<double> >& mat, uint val);
map<uint, double> pasarAMap(const vector<vector<double> >& mat);
vector<double> pasarAVector(const vector<vector<double> >& mat);
VectorMapMatrix  generarRayos(size_t tamMatriz, bool fijos);
VectorMapMatrix  generarRayos_barrido_H(size_t tamMatriz, size_t cada_cuanto);
vector<double> uniformNoise(const vector<double>& t, double init, double end, double sign);
VectorMapMatrix getTraspuesta(const VectorMapMatrix &W);
double ECM(const vector<double>& original, const vector<double>& reconstruido);
pair<vector<double>,short> EG2(vector<vector<double>> &mat, vector<double> bb);
double operator*(const vector<double>& u, const vector<double>& v);
vector<double> operator*(const vector<vector<double> >& M, const vector<double>& v);
void experimentacion_barrido_H(const string& directorio, uint taman_imags, const vector<unsigned short int>& discretizaciones, const vector<pair<float,float> >& ruidos, const vector<unsigned short int>& espacios_entre_censores);
void experimentacion(char tipo, const vector<string>& archivos, string carpeta_salida, uint taman_imags, const vector<unsigned short int>& discretizaciones, const vector<unsigned short int>& cantidades_de_fuentes, const vector<unsigned short int>& separaciones, const vector<pair<float,float> >& ruidos, uint16_t repeticiones);
void listarDirectorio(const string& directorio,  vector<string>& v);
void escribirVector(string nombreArchivo, vector<double>& vector);
void escribirVectorDeVectores(string nombreArchivo, vector<vector<double>>& vector);
void escribirCSV(string nombreArchivo, vector<double>& vector, size_t ancho);

#endif //TC_UTIL_H

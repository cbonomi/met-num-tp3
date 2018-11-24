#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <random>
#include <functional>
#include <tuple>
#include <stdlib.h>
#include "calcular_rayos.h"
#include "util.h"
#include <cmath>
#include <chrono>
#include "VectorMapMatrix.h"


using namespace std;


/**
 * esta funcion toma como parametros las matrices D y V
 * @return el tiempo que tarda la senial en atravesar el cuerpo
 */
vector<double> multMatPorVect(const vector<vector<double> > &M, const vector<double> &v){//recordar que el vector v son las inversas de las velocidades
    const size_t& n = v.size();
    vector<double> res(n);
    for(size_t i = 0; i < n; ++i) {
        res[i] = 0;
        for (size_t k = 0; k < n; ++k)
            res[i] += M[i][k]*v[k];
    }
    return res;
}

vector<vector<double> > trasponer(const vector<vector<double> >& mat){
    const unsigned long& n = mat.size();
    const unsigned long& m = mat[0].size();
    vector<vector<double> > res (m, vector<double>(n));
    for (uint i = 0; i<n;i++)
        for (uint j = 0; j<m;j++)
            res[j][i] = mat[i][j];
    return res;

}

vector<vector<double> > multMatPorMat(const vector<vector<double> >& mat1, const vector<vector<double> >& mat2) {
    const unsigned long& n = mat1.size();
    const unsigned long& m = mat2[0].size();
    const unsigned long& l = mat2.size();
    vector<vector<double> > res (n, vector<double>(m, 0));
    for (uint i = 0; i < n; i++)
        for (uint j = 0; j < m; j++)
            for (uint k = 0; k < l; k++)
                res[i][j] += mat1[i][k]*mat2[k][j];
    return res;
}

VectorMapMatrix multMatPorMat(VectorMapMatrix &mat1, VectorMapMatrix &mat2) {
    	const unsigned long& n = mat1.cantFilas();
    	const unsigned long& m = mat2.cantColumnas();
    	const unsigned long& l = mat2.cantFilas();
    	VectorMapMatrix res = VectorMapMatrix(n, m);
	
    	for (uint i = 0; i < n; i++){
        	for (uint j = i; j < m; j++){
			double acum = 0;
			for (uint k = 0; k < l; k++){
                		acum += mat1.at(i,k)*mat2.at(k,j);
			}
			res.asignar(i,j,acum);
			res.asignar(j,i,acum);
		}
	}
    	return res;
}


vector<vector<uint16_t>> datosAMatriz(uchar &datos, uint ancho, uint alto) {
	vector<vector<uint16_t>> ret (0);
	for (size_t i = 0; i < alto; ++i) {
		vector<uint16_t> fila (0);
		for (size_t j = 0; j < ancho; ++j) {
			int indice = (i*ancho)+j;
			int16_t bt = static_cast<int16_t>(*((char *)&datos + indice));
			fila.push_back(bt);
		}
		ret.push_back(fila);
	}
	return ret;
}

/**
 * Genera Matriz con todos los D_kij (cada fila es una de las matrices D_k).
 * @param tamMatriz tamaño de la imagen discretizada.
 * @param metodo_usado es un numero QUE DEBE VALER 0,1 o 2, y que indica, si es 0, que se usara el metodo de rotaciones
 * iniciando con rayos horizontales, si vale 1, serán unos rayos fijos, que son colocados en los lados horizontales de
 * la imagen y rotaran, si vale 2, estos rayos son colocados en el tope y fondo verticales de la imagen, y tambien rotan,
 * si vale 3 entonces se usa el metodo de horizontales agregando rayos que vengan del tope, si vale 4 usa un metodo en
 * el que evita repetir rayos.
 * @param cantLaseres es la cantidad de laseres que se desean, DEBE SER DIVISOR DE tamMatriz o la función puede tener
 * resultados indeseables, (como minimo puede pasar que no se obtenga la cantidad deseada de laseres, o cosas peores).
 * @param saltear_hasta_n es la cantidad de pixeles rotados que saltearemos despues de cada rayo disparado, el minimo
 * valor permitido es 1 (CON 0 SE ROMPE) y aumentar el valor reduce el tiempo de computo, pero tambien reduce la
 * precision.
 * @return La matriz D con todos los D_k.
 */
/*VectorMapMatrix  generarRayos(size_t tamMatriz, int metodo_usado, int cantLaseres, int saltear_hasta_n);*/
/*
vector<vector<double>> obtenerTrayectorias() {

	for (vector<string>::iterator it = lecturas.begin(); it != lecturas.end(); ++it) {

	}
}*/

bool esTraspuesta(VectorMapMatrix &D, VectorMapMatrix &Dt) {
	bool ret = true;
	for (uint i=0; i<D.cantFilas(); i++) {
		for (uint j=0; j<D.cantColumnas(); j++) {
			ret = ret && D[i][j] == Dt[j][i];
		}
	}
	return ret;
}

vector<double> reconstruirCuerpo(string nombreAchivoEntrada, vector<double>* V, uint tamanoDiscretizacion, double inicioRuido, double finRuido, double signoRuido, size_t* ancho) {
	vector<vector<double> >* cuerpo;
	// 1) tomamos la imagen
	cuerpo = leerCSV(nombreAchivoEntrada);

	// 2) la discretizamos
	vector<vector<double> > cuerpoDiscretizado = discretizar(*cuerpo, tamanoDiscretizacion);
	size_t tamMatriz = cuerpoDiscretizado.size();
    *ancho = cuerpoDiscretizado.size();
	// 3) obtenemos D (la matriz con las trayectorias de los rayos
	VectorMapMatrix  D = generarRayos(tamMatriz, 2, tamMatriz, 1); //tamaño discretizado, metodo a utilizar, cantidad de rayos, pixeles salteados-1.
	// 4) pasamos la imagen discretizada a vector
	vector<double> Vtemp = pasarAVector(cuerpoDiscretizado);
	V = &Vtemp;
	// 5) invertimos el vector V
	/*vector<double> Vinv (V->size(), 0);
	for (uint i = 0; i< V->size(); i++){
		if ((*V)[i] != 0){
			Vinv[i] = 1/(*V)[i];
		}
	} No hay que invertir.*/
	// 6) multiplicamos la matriz D por el vector V
	vector<double> T = D*Vtemp;
	// 7) le aplicamos ruido al vector T
	vector<double> Tr = uniformNoise(T, inicioRuido, finRuido, signoRuido);
	// 8) generamos DtD
	VectorMapMatrix Dt = getTraspuesta(D);
	vector<vector<double>> DtD = Dt*D;//multMatPorMat(Dt,D);
	
	// 9) generamos el vector Dt*T
	vector<double> DtT = Dt*Tr;
	// 10) resolvemos el sistema DtDx = DtT con EG
	pair<vector<double>,short> solucion = EG2(DtD, DtT);
	/*vector<double> Check (V->size(), 0);
	for (uint i = 0; i< V->size(); i++){
		if (abs(solucion.first[i]) > 0.00001){
			Check[i] = 1/solucion.first[i];
		}
	} No hay que invertir.*/

	//cout << ECM(*V,solucion.first) << endl;
	// invertir los valores de la solucion y volverlo a pasar a matriz para luego convertirlo en una imagen que podamos ver
	return solucion.first;
}


//------------------------ Parseo de la entrada -------------------------------//

bool contiene(char *argv[], const string *cadena) {
    string param1 = argv[1], param2 = argv[3], param3 = argv[5];//, param4 = argv[5];
    return param1.compare(*cadena) || param2.compare(*cadena) || param3.compare(*cadena); //|| param4.compare(*cadena);
}

string obtener(char *argv[], const string *cadena) {
    string ret;
    string param1 = argv[1], param2 = argv[3], param3 = argv[5];//, param4 = argv[7];

    if (param1.compare(*cadena) == 0) ret = argv[2];
    if (param2.compare(*cadena) == 0) ret = argv[4];
    if (param3.compare(*cadena) == 0) ret = argv[6];
//    if (param4.compare(*cadena) == 0) ret = argv[8];
    return ret;
}

bool obtenerParametros(int argc, char * argv[], string *ruido, string *nombreArchivoEntrada, string *nombreArchivoSalida) {
    bool ret = false;
    const string param1 = "-r", param2 = "-i", param3 = "-o";

    if (argc == 7 && contiene(argv, &param1) && contiene(argv, &param2) && contiene(argv, &param3)) {
        *ruido = obtener(argv, &param1);
        *nombreArchivoEntrada = obtener(argv, &param2);
        *nombreArchivoSalida = obtener(argv, &param3);
        ret = (ruido != NULL && nombreArchivoEntrada != NULL && nombreArchivoSalida != NULL);
    }
    return ret;
}

//------------------------ Parseo de la entrada -------------------------------//

int main(int argc, char * argv[]) {


    string nombreAchivoEntrada;
    string nombreAchivoSalida;
    vector<double>* V;
    uint discretizacion;
    string ruido;
    size_t ancho;

    if (!obtenerParametros(argc, argv, &ruido, &nombreAchivoEntrada, &nombreAchivoSalida)) {
        cout << "Modo de uso: tp3 -r <nivel_ruido> -i <nombre_archivo_entrada> -o <nombre_archivo_salida>\n";
    } else {
        double nivelRuido = atof(ruido.c_str());
        vector<double> reconstruccion = reconstruirCuerpo(nombreAchivoEntrada, V, 32, nivelRuido, nivelRuido, 0, &ancho);

        escribirCSV(nombreAchivoSalida, reconstruccion, ancho);
    }

    return 0;
}

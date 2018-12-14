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

vector<double> reconstruirCuerpo(string nombreAchivoEntrada, vector<double>* V, uint tamanoDiscretizacion,
                            double nivelRuido, string tipoRuido, uint metodo, uint cantidadLasers,
                            uint separacionRayos, size_t* ancho, const double tolerance, const bool debug) {
	vector<vector<double> >* cuerpo;
	// 1) tomamos la imagen
	cuerpo = leerCSV(nombreAchivoEntrada);

	// 2) la discretizamos
	uint granularidad = cuerpo->size()/tamanoDiscretizacion; 
	vector<vector<double> > cuerpoDiscretizado = discretizar(*cuerpo, granularidad);
	size_t tamMatriz = cuerpoDiscretizado.size();
    *ancho = cuerpoDiscretizado.size();
	// 3) obtenemos D (la matriz con las trayectorias de los rayos
	VectorMapMatrix  D = generarRayos(tamMatriz, metodo, cantidadLasers, separacionRayos); //tamaño discretizado, metodo a utilizar, cantidad de rayos, pixeles salteados-1.
	// 4) pasamos la imagen discretizada a vector
	vector<double> Vtemp = pasarAVector(cuerpoDiscretizado);
	//V = &Vtemp;
	
	// 6) multiplicamos la matriz D por el vector V
	vector<double> T = D*Vtemp;
	// 7) le aplicamos ruido al vector T
    //copiamos T a Tr
    vector<double> Tr;
    for (uint i=0;i<T.size();i++)
        Tr.push_back(T[i]);

    if (tipoRuido.compare("A") == 0) // aplicamos ruido aditivo
	    Tr = AWGNNoise(T, Vtemp, nivelRuido);
	if (tipoRuido.compare("M") == 0) // aplicamos ruido multiplicativo
        Tr = MWGNNoise(T, Vtemp, nivelRuido);

	// 8) generamos DtD
	VectorMapMatrix Dt = getTraspuesta(D);
	vector<vector<double>> DtD = Dt*D;//multMatPorMat(Dt,D);
	
	// 9) generamos el vector Dt*T
	//vector<double> DtT = Dt*Tr;

	// 10) resolvemos el sistema DtDx = DtT con CML
	vector<vector<double> > A = D.convert_to_vec_matrix();
	vector<double> solucion = CML(A, Tr,tolerance,debug);

    return solucion;
}


//------------------------ Parseo de la entrada -------------------------------//

bool contiene(char *argv[], const string *cadena) {
    string param1 = argv[1], param2 = argv[3], param3 = argv[5], param4 = argv[7], param5 = argv[9];
    return param1.compare(*cadena) || param2.compare(*cadena) || param3.compare(*cadena) || param4.compare(*cadena)
            || param5.compare(*cadena);
}

string obtener(char *argv[], const string *cadena) {
    string ret;
    string param1 = argv[1], param2 = argv[3], param3 = argv[5], param4 = argv[7], param5 = argv[9];

    if (param1.compare(*cadena) == 0) ret = argv[2];
    if (param2.compare(*cadena) == 0) ret = argv[4];
    if (param3.compare(*cadena) == 0) ret = argv[6];
    if (param4.compare(*cadena) == 0) ret = argv[8];
    if (param5.compare(*cadena) == 0) ret = argv[10];

    return ret;
}

bool obtenerParametros(int argc, char * argv[], string *ruido, string *nombreArchivoEntrada,
                        string *nombreArchivoSalida, string *tipoRuido, string *discretizacion) {
    bool ret = false;
    const string param1 = "-r", param2 = "-i", param3 = "-o", param4 = "-t", param5 = "-d";

    if (argc == 11 && contiene(argv, &param1) && contiene(argv, &param2) && contiene(argv, &param3)
        && contiene(argv, &param4) && contiene(argv, &param5)) {
        *ruido = obtener(argv, &param1);
        *nombreArchivoEntrada = obtener(argv, &param2);
        *nombreArchivoSalida = obtener(argv, &param3);
        *tipoRuido = obtener(argv, &param4);
        *discretizacion = obtener(argv, &param5);
        ret = (ruido != NULL && nombreArchivoEntrada != NULL && nombreArchivoSalida != NULL && tipoRuido != NULL
                && discretizacion != NULL);
    }
    return ret;
}

//------------------------ Parseo de la entrada -------------------------------//

int main(int argc, char * argv[]) {

    string nombreArchivoEntrada;
    string nombreArchivoSalida;
    vector<double>* V;
    uint discretizacion = 16;
    double tolerance = 0.0;
    double nivelRuido = 0.0;
    string tipoRuido = "A";
    size_t ancho;
    uint metodo = 2;
    uint separacionRayos = 1;
    uint cantidadRayos = 0;

    bool debug =false;
    bool in = false;
    bool out = false;
    for (int c=1;c<argc;++c){
        string arg = argv[c];
        if (arg=="-o"){
            out = true;
            nombreArchivoSalida = argv[c+1];
        }
        if (arg=="-i"){
            in = true;
            nombreArchivoEntrada = argv[c+1];
        }
        if (arg=="-r"){
            nivelRuido = atof(argv[c+1]);
        }
        if (arg=="-t"){
            tipoRuido = argv[c+1];
        }
        if (arg=="-d"){
            discretizacion = atoi(argv[c+1]);
        }
        if (arg=="-nt"){
            tolerance = atof(argv[c+1]);
        }
        if (arg=="-debug"){
            debug = true;
        }
        /*if (arg=="-numcond" or arg=="-nc"){
            debug = true;
        }*/
        if (arg=="-method" or arg=="-m"){
            metodo = atoi(argv[c+1]);
        }
        if (arg=="-s"){
            separacionRayos = atoi(argv[c+1]);
        }
        if (arg=="-cantRayos" or arg=="-cr"){
            cantidadRayos = atoi(argv[c+1]);
        }
    }

    if (cantidadRayos==0){
        //Valor por defecto
        cantidadRayos = discretizacion;
    }
    
    bool faltan_parametros = not (in and out);
    if (faltan_parametros){
        cout << "Modo de uso: tp3 -r <nivel_ruido> -t <tipo_ruido> -i <nombre_archivo_entrada> -o <nombre_archivo_salida>\n";
        cout << "Las opciones para tipo de ruido son:" << endl;
        cout << "\t M: ruido multiplicativo" << endl;
        cout << "\t A: ruido aditivo." << endl;
        cout << "---- Parametros opcionales: ----" << endl;
        cout << "-d <tam_discretizacion> : 16 por defecto." << endl;
        cout << "-nt <nivel_tolerancia> : Valor por defecto 0" << endl;
        cout << "     t=0 : No se filtran valores singulares" << endl;
        cout << "     t=0.001 Se filtran valores singulares 1000 veces mas chicos que el valor singular de mayor modulo" << endl;
        //cout << "-nc/-numcond : Imprime el numero de condicion de la matriz" << endl;
        cout << "-m <num_metodo> selecciona el metodo de creacion de rayos, 1 por defecto" << endl;
        cout << "-s <num_separacion> selecciona el espaciado entre cada rayo generado, 1 por defecto (sin espaciado)" << endl;
        cout << "-l <num_rayos> selecciona la cantidad de rayos, por defecto el tamano completo" << endl;
        
    } else {
        vector<double> reconstruccion = reconstruirCuerpo(nombreArchivoEntrada, V, discretizacion, nivelRuido,tipoRuido, 
                                        metodo, cantidadRayos, separacionRayos, &ancho,tolerance,debug);
        escribirCSV(nombreArchivoSalida, reconstruccion, ancho);
    }

    return 0;
}

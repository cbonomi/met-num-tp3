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

typedef unsigned char uchar;

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

void imprimirVector(const vector<unsigned short> & vect, uint tam){
    for(uint i = 0; i< tam; i++){
        for(uint j = 0; j < tam; j++){
            cout << vect[i*tam+j] << ",";
        } cout << endl;
    }

}

double reconstruirCuerpo(string nombreAchivoEntrada, string nombreAchivoSalida, vector<double>* V, uint metodo, uint tamanoDiscretizacion, double nivelRuido, size_t* ancho, int cantidadLaseres, int separacionRayos) {
    vector<vector<double> >* cuerpo;
    // 1) tomamos la imagen
    cuerpo = leerCSV(nombreAchivoEntrada);

    // 2) la discretizamos
    vector<vector<double> > cuerpoDiscretizado = discretizar(*cuerpo, tamanoDiscretizacion);
    size_t tamMatriz = cuerpoDiscretizado.size();
    *ancho = cuerpoDiscretizado.size();
    // 3) obtenemos D (la matriz con las trayectorias de los rayos
    VectorMapMatrix D = generarRayos(tamMatriz, metodo, cantidadLaseres, separacionRayos); //tamaño discretizado, metodo a utilizar, cantidad de rayos, pixeles salteados-1.
    // 4) pasamos la imagen discretizada a vector
    vector<double> Vtemp = pasarAVector(cuerpoDiscretizado);
    V = &Vtemp;

    // 6) multiplicamos la matriz D por el vector V
    vector<double> T = D*Vtemp;
    // 7) le aplicamos ruido al vector T
    vector<double> Tr = MWGNNoise(T, nivelRuido);
/*
    // 8) generamos DtD
    VectorMapMatrix Dt = getTraspuesta(D);
    vector<vector<double>> DtD = Dt*D;//multMatPorMat(Dt,D);

    // 9) generamos el vector Dt*T
    //vector<double> DtT = Dt*Tr;

    // 10) resolvemos el sistema DtDx = DtT con CML
    vector<vector<double> > A = D.convert_to_vec_matrix();
    vector<double> solucion = CML(A, Tr);

    // 10) resolvemos el sistema DtDx = DtT con EG
    //pair<vector<double>,short> solucion = EG2(DtD, DtT);
*/


    // 8) generamos DtD
    VectorMapMatrix Dt = getTraspuesta(D);
    vector<vector<double>> DtD = Dt*D;//multMatPorMat(Dt,D);

    // 9) generamos el vector Dt*T
    vector<double> DtT = Dt*Tr;

    // 10) resolvemos el sistema DtDx = DtT con EG
    pair<vector<double>,short> solucion = EG2(DtD, DtT);

    //cout << ECM(*V,solucion.first) << endl;
    // invertir los valores de la solucion y volverlo a pasar a matriz para luego convertirlo en una imagen que podamos ver
    //string salida = nombreAchivoEntrada + "disc.csv";
    escribirCSV(nombreAchivoSalida.c_str(), solucion.first, *ancho);
//    cout << nombreAchivoSalida << endl;
    return calcularPSNR(Vtemp, solucion.first);
}


//------------------------ Parseo de la entrada -------------------------------//

bool contiene(char *argv[], const string *cadena) {
    string param1 = argv[1], param2 = argv[3], param3 = argv[5], param4 = argv[7], param5 = argv[9], param6 = argv[11],
            param7 = argv[13];
    return param1.compare(*cadena) || param2.compare(*cadena) || param3.compare(*cadena) || param4.compare(*cadena)
           || param5.compare(*cadena) || param6.compare(*cadena) || param7.compare(*cadena);
}

string obtener(char *argv[], const string *cadena) {
    string ret;
    string param1 = argv[1], param2 = argv[3], param3 = argv[5], param4 = argv[7], param5 = argv[9], param6 = argv[11],
            param7 = argv[13];

    if (param1.compare(*cadena) == 0) ret = argv[2];
    if (param2.compare(*cadena) == 0) ret = argv[4];
    if (param3.compare(*cadena) == 0) ret = argv[6];
    if (param4.compare(*cadena) == 0) ret = argv[8];
    if (param5.compare(*cadena) == 0) ret = argv[10];
    if (param6.compare(*cadena) == 0) ret = argv[12];
    if (param7.compare(*cadena) == 0) ret = argv[14];
    return ret;
}

bool obtenerParametros(int argc, char * argv[], string *ruido, string *nombreArchivoEntrada,
                       string *nombreArchivoSalida, string *metodo, string *discretizacion,
                       string *cantidadLasers, string *separacionRayos) {
    bool ret = false;
    const string param1 = "-r", param2 = "-i", param3 = "-o", param4 = "-m", param5 = "-d", param6 = "-l", param7 = "-s";

    if (argc == 15 && contiene(argv, &param1) && contiene(argv, &param2) && contiene(argv, &param3) &&
        contiene(argv, &param4) && contiene(argv, &param5) && contiene(argv, &param6) && contiene(argv, &param7)) {
        *ruido = obtener(argv, &param1);
        *nombreArchivoEntrada = obtener(argv, &param2);
        *nombreArchivoSalida = obtener(argv, &param3);
        *metodo = obtener(argv, &param4);
        *discretizacion = obtener(argv, &param5);
        *cantidadLasers = obtener(argv, &param6);
        *separacionRayos = obtener(argv, &param7);
        ret = (ruido != NULL && nombreArchivoEntrada != NULL && nombreArchivoSalida != NULL && metodo != NULL
               && discretizacion != NULL && cantidadLasers != NULL && separacionRayos != NULL);
    }
    return ret;
}

//------------------------ Parseo de la entrada -------------------------------//

int main(int argc, char * argv[]) {


    string nombreAchivoEntrada;
    string nombreAchivoSalida;
    vector<double>* V;
    string ruido;
    size_t ancho;
    string m;
    string d;
    string l;
    string s;

    if (!obtenerParametros(argc, argv, &ruido, &nombreAchivoEntrada, &nombreAchivoSalida, &m, &d, &l, &s)) {
        cout << "Modo de uso: tp3 -r <nivel_ruido> -i <nombre_archivo_entrada> -o <nombre_archivo_salida> -m <metodo> "
                "-d <discretizacion> -l <cantidad_lasers> -s <separacion>\n";
    } else {
        double nivelRuido = atof(ruido.c_str());
        int metodo = atoi(m.c_str());
        uint discretizacion = atoi(d.c_str());
        int cantidadLasers = atoi(l.c_str());
        int separacionRayos = atoi(s.c_str());
//        cout << nombreAchivoEntrada << endl << nombreAchivoSalida << endl;
        double PSNR = reconstruirCuerpo(nombreAchivoEntrada, nombreAchivoSalida, V, metodo, discretizacion, nivelRuido, &ancho,
                                        cantidadLasers, separacionRayos);
        cout << PSNR;
//        escribirCSV(nombreAchivoSalida, reconstruccion, ancho);
    }

    return 0;
}


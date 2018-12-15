#include "util.h"

using namespace std;

const vector<string> explode(const string& s, const char& c)
{
    string buff{""};
    vector<string> v;

    for(auto n:s)
    {
        if(n != c) buff+=n; else
        if(n == c && buff != "") { v.push_back(buff); buff = ""; }
    }
    if(buff != "") v.push_back(buff);

    return v;
}


vector<double> convertirAVectorDeValores(vector<string> lecturas) {
    vector<double>* ret = new vector<double>(0);
    for (vector<string>::iterator it = lecturas.begin() ; it != lecturas.end(); ++it) {
        string lectura = *it;
        ret->push_back(stoi(lectura));
    }
    return *ret;
}

vector<vector<double>>* leerCSV(string nombreArchivo) {
    fstream entrada(nombreArchivo, ios_base::in);
    if(entrada.fail()){
        throw runtime_error("no se existe el archivo " + nombreArchivo + "!");
    }
    vector<vector<double>>* ret = new vector<vector<double>>(0);

    string lectura;
    vector<string> lecturas;
    bool path = true;
    while(entrada >> lectura) {
        ret->push_back(convertirAVectorDeValores(explode(lectura, ',')));
    }
    entrada.close();
    return ret;
}




vector<vector<double>> discretizar(const vector<vector<double> >& mat, uint val){//supongo que la matriz mat es cuadrada
	vector<vector<double>> res (mat.size()/val, vector<double> (mat.size()/val));
	for(uint i = 0; i< res.size(); i++){
		for(uint j = 0; j < res.size(); j++){
			double temp = 0;
			for (uint k = 0+i*val; k < (i+1)*val; k++){
				for(uint h = 0+j*val; h < (j+1)*val; h++){
					temp+= mat[k][h];
				}
			}
			res[i][j] = temp/(val*val);
		}
	}
	return res;
}

map<uint, double> pasarAMap(const vector<vector<double> >& mat){
    map<uint, double> res;
    for(uint i = 0; i< mat.size(); i++){
        for(uint j = 0; j < mat[0].size(); j++){
            if(abs(mat[i][j]) > 0.001) {
                res[i*mat[0].size()+j] = mat[i][j];
            }
        }
    }
    return res;

}

vector<double> pasarAVector(const vector<vector<double> >& mat){
	vector<double> res (mat.size()*mat[0].size(),0);
	for(uint i = 0; i< mat.size(); i++){
		for(uint j = 0; j < mat[0].size(); j++){
			res[i*mat[0].size()+j] = mat[i][j];
		}
	}
	return res;

}

VectorMapMatrix generarRayos_barrido_H(size_t tamMatriz, size_t cada_cuanto) {
    vector<pair<uint,uint> > laseres = crearLaseres(tamMatriz, cada_cuanto, cada_cuanto/2, 0);
    vector<pair<uint,uint> > sensores = crearPuntosDeFin(laseres, tamMatriz);
    VectorMapMatrix D_ks(0, tamMatriz*tamMatriz);
    while(sensores[0].first != tamMatriz - 1 or sensores[0].second != 0) { //Esto quiza es dificil de ver, pero para los laseres izquierdos
        // que se saltean el horizontal, este es la ultima posicion interesante a la que apuntan. NOTA IMPORTANTE,
        // SI SE HACEN MAS DE UN SALTO PUEDE QUE ESTO NO TERMINE. ASIQUE CUIDADO CON PONER MAS DE UN rotarLaseres.
        for(uint i = 0; i < laseres.size(); i++) {
            map<uint, double> D_k_map = pasarAMap(trazar_recta_en_matriz_D(laseres[i], sensores[i], tamMatriz));
            D_ks.agregarFila(D_k_map);
        }
        barrerLaseres_H(laseres,sensores,tamMatriz);
    }
    return D_ks;
}

vector<double> uniformNoise(const vector<double>& t, double init, double end, double sign){
    vector<double> res(t.size());
    default_random_engine generator;
    uniform_real_distribution<double> distribution(init,end);
    for(uint i = 0; i< t.size(); i++){
        double number = distribution(generator);
        if (sign != 0){
            res[i] = sign*number*t[i] + t[i];
        }
        else {
            random_device rd;
            mt19937 gen(rd());
            uniform_int_distribution<> dis(1, 2);
            if (dis(gen) == 1){
                res[i] = number*t[i] + t[i];
            }
            else {
                res[i] = t[i] - number*t[i];
            }
        }
    }
    return res;
}

VectorMapMatrix getTraspuesta(const VectorMapMatrix &W) {
    VectorMapMatrix ret(W.cantColumnas(), W.cantFilas());

    for(uint i = 0; i < W.cantFilas(); ++i)
        for (unsigned int j=0; j<W.cantColumnas(); ++j)
            ret.asignar(j, i, W.at(i, j));
    return ret;

}

vector<double> CML(vector<vector<double>> &mat, vector<double> b,const double tolerance,const bool debug) {
    vector<vector<double> > Ut;
    vector<double> Sigma;
    vector<vector<double> > Vt;
    calcular_svd(mat,Ut,Sigma,Vt,tolerance,debug);
    vector<double> b_prima = mult_matr_por_vect(Ut,b);
    vector<double> res = vector<double>(mat.size(),0.0);
    
    for (uint i = 0; i < Sigma.size(); ++i){
       // TODO: formula:
        vector<double> vec = Vt[i];
        multVecEsc(vec,b_prima[i]);
        multVecEsc(vec,(1/Sigma[i]));

        res = sumaVec(res,vec);

        //res = res + (b_prima[i] / Sigma[i]) * Vt.fila(i);
    }
    return res;
}




double operator*(const vector<double>& u, const vector<double>& v){   //Deben ser del mismo tamaño.
    double res = 0;
    for(size_t i = 0; i < u.size(); ++i)
        res += u[i]*v[i];
    return res;
}

vector<double> operator*(const vector<vector<double> >& M, const vector<double>& v){
    vector<double> res(v.size());
    for(size_t i = 0; i < M.size(); ++i)
        res[i] = M[i]*v;
    return res;
}

void listarDirectorio(const string& directorio,  vector<string>& v)
{
    string nomArch;
    DIR* dirp = opendir(directorio.c_str());
    struct dirent * dp;
    if (dirp == NULL) {
        throw runtime_error("no se encontro directorio " + directorio + "!");
    }
    while ((dp = readdir(dirp)) != NULL) {
        string nomArch = dp->d_name;
        if (nomArch.compare(".") != 0 && nomArch.compare("..") != 0)
            v.push_back(directorio + "/" + nomArch);
    }
    closedir(dirp);
}

void escribirVector(string nombreArchivo, vector<double>& vector) {
    ofstream salida(nombreArchivo, ios_base::out);
    for (uint i=0; i < vector.size(); i++) {
        salida << vector[i] << endl;
    }
    salida.close();
}


void escribirVectorDeVectores(string nombreArchivo, vector<vector<double>>& vector) {
    ofstream salida(nombreArchivo, ios_base::out);
    string linea = "";
    for (uint i = 0; i < vector.size(); i++) {
        for (uint j = 0; j < vector[i].size(); j++) {
            linea += to_string(vector[i][j]) + " ";

        }
        salida << linea << endl;
        linea = "";
    }
    salida.close();
}

void escribirCSV(string nombreArchivo, vector<double>& vector, size_t ancho) {
    ofstream salida(nombreArchivo, ios_base::out);
    string linea = "";
    double valor;
    for (uint j=0; j<ancho; j++) {
        for (uint i = 0; i < ancho-1; i++) {
            valor = floor(vector[i + j*ancho]);
            linea += to_string((signed short) valor) + ",";
        }
        valor = floor(vector[ancho-1 + j*ancho]);
        linea += to_string((signed short) valor);
        salida << linea << endl;
        linea = "";
    }
    salida.close();
}

double WGN_generate()
{/* Generates additive white Gaussian Noise samples with zero mean and a standard deviation of 1. */

    double temp1;
    double temp2;
    double result;
    int p;

    p = 1;

    while( p > 0 )
    {
        temp2 = ( rand() / ( (double)RAND_MAX ) ); /*  rand() function generates an
                                                       integer between 0 and  RAND_MAX,
                                                       which is defined in stdlib.h.
                                                   */

        if ( temp2 == 0 )
        {// temp2 is >= (RAND_MAX / 2)
            p = 1;
        }// end if
        else
        {// temp2 is < (RAND_MAX / 2)
            p = -1;
        }// end else

    }// end while()

    temp1 = cos( ( 2.0 * (double)PI ) * rand() / ( (double)RAND_MAX ) );
    result = sqrt( -2.0 * log( temp2 ) ) * temp1;

    return result;	// return the generated random sample to the caller

}// end AWGN_generator()


double calcularMedia(const vector<double>& t) {
    double res = 0;
    for (const auto &valor : t) {
        res += valor;
    }
    return res/t.size();
}

double calcularDesvio(const vector<double>& t){
    double res(t.size());

    double media = calcularMedia(t);

    for (const auto &valor : t) {
        res += pow(valor - media, 2);
    }
    return sqrt(res/t.size()-1);
}

/*
 * Devuelve un vector de tamaño n con el ruido (basado en el desvio pasado como parametro) generado
 */
vector<double> WGNNoise(size_t n, double desvio){
    vector<double> res(n);
    for(uint i = 0; i< n; i++){
        double noise = WGN_generate()*desvio;
        res[i] = noise;
    }
    return res;
}


/*
 * Dado un vector calcula su desvio y genera otro agregandole ruido multiplicativo
 * de acuerdo al porcentaje del desvio standar del vector pasado como parametro
 */
vector<double> MWGNNoise(const vector<double>& t, const vector<double>& imagen, double porcentajeDeRuido){
    uint n = t.size();
    vector<double> res(n);
    double desvio = calcularDesvio(imagen) * porcentajeDeRuido;
    for(uint i = 0; i< n; i++){
        double noise = WGN_generate()*desvio;
        res[i] = t[i] * noise;
    }
    return res;
}

/*
 * Dado un vector calcula su desvio y genera otro agregandole ruido aditivo
 * de acuerdo al porcentaje del desvio standar del vector pasado como parametro
 */
vector<double> AWGNNoise(const vector<double>& t, const vector<double>& imagen, double porcentajeDeRuido){
    uint n = t.size();
    vector<double> res(n);
    double desvio = calcularDesvio(imagen) * porcentajeDeRuido;
    for(uint i = 0; i< n; i++){
        double noise = WGN_generate()*desvio;
        res[i] = t[i] + noise;
    }
    return res;
}

double ECM(const vector<double>& original, const vector<double>& reconstruido) {
    uint n = original.size();
    double ret = 0;
    double dif;
    int val;
    for(uint i = 0; i< n; i++){
        if (reconstruido[i] < 0)
            val = 0;
        else
        if (reconstruido[i] > 255)
            val = 255;
        else
            val = reconstruido[i];
        dif = original[i] - val;
        ret += dif*dif;
    }
    return ret/n;
}

long double calcularPSNR(const vector<double>& original, const vector<double>& reconstruido) {
    return 10 * log10 (MAX_u_cuadrado/ECM(original, reconstruido));
}

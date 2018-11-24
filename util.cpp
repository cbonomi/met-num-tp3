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

double ECM(const vector<double>& original, const vector<double>& reconstruido) {
    uint n = original.size();
    double ret = 0;
    double dif;
    for(uint i = 0; i< n; i++){
        dif = original[i] - reconstruido[i];
        ret += dif*dif;
    }
    return ret/(n*n);
}

pair<vector<double>,short> EG2(vector<vector<double>> &mat, vector<double> bb) {
	unsigned int i,j,l;
	vector<double> res(mat[0].size(),0);
	short status = 0; //status default, el sistema tiene una unica solucion posible
	double A_kk, A_jk;
	vector<vector<double>> copy = mat;
	vector<double> copy2 = bb;
	bool cont;


	for(i = 0; i < copy[0].size()-1; i++){ //itero sobre las filas, excepto la ultima porque ahi no tengo que hacer nada
		cont = false;
		for(j = i; j < copy.size(); j++){ //itero sobre las filas desde i en adelante, estaria por fijarme si tengo que hacer o no calculo en el paso i de la EG
			if(abs(copy[j][i]) > 0.00001){ //si no hay un 0 en la posicion j,i
				cont = true;
				if(abs(copy[i][i]) <= 0.00001){
					copy[i].swap(copy[j]); //cambio de lugar las filas porque habia un 0 en la diagonal pero no en el resto de la columna
                    			double temp = copy2[i];
                   			copy2[i] = copy2[j];         //como se cambiaron de lugar las filas, también se cambian de lugar los valores de "bb"
                    			copy2[j] = temp;
                		}
				break;
			}
		}
		A_kk = copy[i][i];
		for(j = i+1; j < mat.size(); j++){ //cálculo del paso i si corresponde
			
			if (!cont){break;} //si me tengo que saltear este paso no calculo nada
			if(abs(copy[j][i]) > 0.00001){//si el elemento j,i es 0 no hago nada en la fila j
				A_jk = copy[j][i];
				double temp3 = -A_jk/A_kk;
				if (abs(temp3) > 0.00001){
				for(l = i+1; l < mat[0].size(); l++){
					double temp = temp3*copy[i][l];
					if (abs(temp) > 0.00001 and abs(copy[i][l]) > 0.00001) 
					copy[j][l] +=temp;
					if (abs(copy[j][l]) <= 0.00001)
					copy[j][l] = 0;
				}
				double temp2 = copy2[i]*temp3;
				if (abs(temp2) > 0.00001 and abs(copy2[i]) > 0.00001)
				copy2[j] += temp2;
				if (abs(copy2[j]) <= 0.00001)
				copy2[j] = 0;
				}
			} //no me olvido de actualizar el vector b
		}
		
	}
	
	


	for(i = 0; i < copy[0].size(); i++){
		j = copy[0].size()-1-i;
		if(copy[j][j] == 0 && copy2[j] != 0){
			status = -1; //el sistema es incompatible
			break;
		}
		if(copy[j][j] == 0 && copy2[j] == 0){
			status = 1; //hay infinitos resultados
			res[j] = 0;
		}
		else{
			res[j] = copy2[j]/copy[j][j]; //tengo A_jj*x_j = b_j, paso dividiendo el A_jj
			
			if (j!=0){
				for(unsigned int l = 0; l < j; l++){
					if (abs(copy[l][j]) > 0.00001)
					copy2[l] -= res[j]*copy[l][j]; //esto es importante, al b_l con l de 0 a j-1 le paso restando el A_lj*x_j, porque ya conozco el resultado de X_j, de forma que en la siguiente iteracion solo voy a tener algo de esta pinta A_jj*x_j = b_j
				}
			}
		}
	}
	return make_pair(res,status);
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

/**
 * @param directorio: nombre del directorio con las imágenes que se usaran.
 * @param taman_imags: cantidad de pixeles de las imagenes.
 * @param discretizacion: vector con las distintas cantidades de pixeles de lado y alto por casillero.
 * @param ruido: vector con los distintos intervalos del porcentaje de ruido (expresado como valor entre 0 y 1)
 * @param espacio_entre_censores
 */
void experimentacion_barrido_H(const string& directorio, uint taman_imags, const vector<unsigned short int>& discretizaciones, const vector<pair<float,float> >& ruidos, const vector<unsigned short int>& espacios_entre_censores) {   //Necesito saber el tamaño de las imagenes de antemano.
    vector<string> archivos;
    listarDirectorio(directorio, archivos);
    /*archivos.push_back("Imagenes_para_probar/1.2.826.0.1.3680043.2.656.1.138.1.csv");*/
    ofstream salida;
    for(size_t ind_disc = 0; ind_disc < discretizaciones.size(); ++ind_disc){
        for(size_t ind_espac = 0; ind_espac < espacios_entre_censores.size(); ++ind_espac){
            if(taman_imags/discretizaciones[ind_disc] > espacios_entre_censores[ind_espac]/2){ //Si hay espacio para al menos 1 fuente de lasers
                VectorMapMatrix D = generarRayos_barrido_H(taman_imags/discretizaciones[ind_disc], espacios_entre_censores[ind_espac]);
                VectorMapMatrix Dt = getTraspuesta(D);
                vector<vector<double> > Dt_D = Dt * D;
                salida.open("resultados de prueba/Discretizacion:"+to_string(discretizaciones[ind_disc])+" espaciado:"+to_string(espacios_entre_censores[ind_espac])+" .txt");
                salida.close(); //La intención de estas 2 lineas es poner en blanco el archivo si ya existe.
                for(size_t ind_arch = 0; ind_arch < archivos.size(); ++ind_arch){
                    vector<vector<double> > *imagen_entera = leerCSV(archivos[ind_arch]);
                    vector<vector<double> > imagen_discreta = discretizar(*imagen_entera, discretizaciones[ind_disc]);
                    vector<double> vec_imagen_discreta = pasarAVector(imagen_discreta);
                    vector<double> t_sin_ruido = D * vec_imagen_discreta;
                    salida.open("resultados de prueba/Discretizacion:"+to_string(discretizaciones[ind_disc])+" espaciado:"+to_string(espacios_entre_censores[ind_espac])+" .txt", ios::app);
                    salida << "Imagen "+archivos[ind_arch]+":\t";
                    salida.close();
                    for(size_t ind_ruido = 0; ind_ruido < ruidos.size(); ++ind_ruido){
                        vector<double> t_con_ruido = uniformNoise(t_sin_ruido, ruidos[ind_ruido].first, ruidos[ind_ruido].second, 0);
                        pair<vector<double>, short> v = EG2(Dt_D, Dt * t_con_ruido);
                        double error = ECM(vec_imagen_discreta, v.first);
                        salida.open("resultados de prueba/Discretizacion:"+to_string(discretizaciones[ind_disc])+" espaciado:"+to_string(espacios_entre_censores[ind_espac])+" .txt", ios::app);
                        salida << error << ",\t";
                        salida.close();
                    }
                    salida.open("resultados de prueba/Discretizacion:"+to_string(discretizaciones[ind_disc])+" espaciado:"+to_string(espacios_entre_censores[ind_espac])+" .txt", ios::app);
                    salida << endl;
                    salida << endl; //Lo hago 2 veces para mejor visibilidad.
                    salida.close();
                    delete imagen_entera;
                }
            }
        }
    }
}

//#define repeticiones 20

/**
 * @param tipo: Si es 'H' se hace un barrido horizontal.
 *              Si es 'V' se hace un barrido vertical.
 *              Si es 'O' (promedio entre 'H' y 'V') se hace un barrido horizontal y vertical.
 *              Si es 'h' se hace un barrido horizontal sin repeticiones.
 *              Si es 'v' se hace un barrido vertical sin repeticiones.
 *              Si es 'o' se hace un barrido horizontal y vertical con menos repeticiones.
 *              Si es 'r' se usan rotaciones.
 * @param archivos: vector con los distintos nombres de las imágenes que se usaran.
 * @param carpeta_salida: nombre de la carpeta donde se guardarán los resultados.
 * @param taman_imags: cantidad de pixeles de las imagenes.
 * @param discretizacion: vector con las distintas cantidades de pixeles de lado y alto por casillero.
 * @param cantidades_de_fuentes: vector con las distintas cantidades de fuentes de rayos.
 * @param separaciones: vector con las distintas separaciones entre rayos.
 * @param ruido: vector con los distintos intervalos del porcentaje de ruido (expresado como valor entre 0 y 1)
 * @param repeticiones: cantidad de veces que se repite cada parte del proceso cuyo tiempo queremos medir.
 */
void experimentacion(char tipo, const vector<string>& archivos, string carpeta_salida, uint taman_imags, const vector<unsigned short int>& discretizaciones, const vector<unsigned short int>& cantidades_de_fuentes, const vector<unsigned short int>& separaciones, const vector<pair<float,float> >& ruidos, uint16_t repeticiones) {   //Necesito saber el tamaño de las imagenes de antemano.
    ofstream salida;
    for(size_t ind_disc = 0; ind_disc < discretizaciones.size(); ++ind_disc){
        for(size_t ind_fuent = 0; ind_fuent < cantidades_de_fuentes.size(); ++ind_fuent){
            for(size_t ind_separ = 0; ind_separ < separaciones.size(); ++ind_separ){
                uint cant_casilleros = taman_imags/discretizaciones[ind_disc];
                if(cantidades_de_fuentes[ind_fuent] <= cant_casilleros){    // && separaciones[ind_separ] < cant_casilleros/2){ //Quiero que cada fuente genere al menos 6 o 4 rayos aproximadamente (6 para blos barridos y 4 para la rotación)
                    VectorMapMatrix D;
                    string nombre_arch_salida;
                    unsigned long comienzo, final, ciclos_clock[repeticiones];
                    for(uint8_t i = 0; i < repeticiones; ++i){
                        if(tipo == 'r') {  //Rotaciones
                            RDTSC_START(comienzo);
                            D = generarRayos(cant_casilleros, 0, cantidades_de_fuentes[ind_fuent], separaciones[ind_separ]);
                            RDTSC_STOP(final);
                            nombre_arch_salida = carpeta_salida + "/Tipo:R";
                        }else if(tipo == 'H'){  //Barrido horizontal
                            RDTSC_START(comienzo);
                            D = generarRayos(cant_casilleros, 1, cantidades_de_fuentes[ind_fuent], separaciones[ind_separ]);
                            RDTSC_STOP(final);
                            nombre_arch_salida = carpeta_salida + "/Tipo:H";
                        }else if(tipo == 'V') {  //Barrido vertical
                            RDTSC_START(comienzo);
                            D = generarRayos(cant_casilleros, 2, cantidades_de_fuentes[ind_fuent], separaciones[ind_separ]);
                            RDTSC_STOP(final);
                            nombre_arch_salida = carpeta_salida + "/Tipo:V";
                        }else if(tipo == 'O') {  //Barrido vertical y horizontal
                            RDTSC_START(comienzo);
                            D = generarRayos(cant_casilleros, 3, cantidades_de_fuentes[ind_fuent], separaciones[ind_separ]);
                            RDTSC_STOP(final);
                            nombre_arch_salida = carpeta_salida + "/Tipo:HyV";
                        }else if(tipo == 'o') {  //Barrido vertical y horizontal con menos repeticiones.
                            RDTSC_START(comienzo);
                            D = generarRayos(cant_casilleros, 4, cantidades_de_fuentes[ind_fuent], separaciones[ind_separ]);
                            RDTSC_STOP(final);
                            nombre_arch_salida = carpeta_salida + "/Tipo:HyV_SR";
                        }else if(tipo == 'h') {  //Barrido horizontal sin repeticiones.
                            RDTSC_START(comienzo);
                            D = generarRayos(cant_casilleros, 5, cantidades_de_fuentes[ind_fuent], separaciones[ind_separ]);
                            RDTSC_STOP(final);
                            nombre_arch_salida = carpeta_salida + "/Tipo:H_SR";
                        }else if(tipo == 'v') {  //Barrido vertical sin repeticiones.
                            RDTSC_START(comienzo);
                            D = generarRayos(cant_casilleros, 6, cantidades_de_fuentes[ind_fuent], separaciones[ind_separ]);
                            RDTSC_STOP(final);
                            nombre_arch_salida = carpeta_salida + "/Tipo:V_SR";
                        }
                        ciclos_clock[i] = final - comienzo;
                    }
                    nombre_arch_salida += " Discretizacion:"+to_string(discretizaciones[ind_disc])+" cantidad_fuentes:"+to_string(cantidades_de_fuentes[ind_fuent])+" separacion:"+to_string(separaciones[ind_separ])+" .txt";
                    salida.open(nombre_arch_salida);
                    salida << "Cantidad rayos: " << D.cantFilas() << endl;
                    salida << "Cantidad ciclos del calculo de los mismos: ";
                    for(uint8_t i = 0; i < repeticiones; ++i){
                        salida << "N°" << i+1 << " " << ciclos_clock[i] << "; ";
                    }
                    salida << endl << endl;
                    salida.close();
                    VectorMapMatrix Dt;
                    vector<vector<double> > Dt_D;
                    unsigned long ciclos_antes_de_leer_imagenes[repeticiones];
                    for(uint8_t i = 0; i < repeticiones; ++i) {
                        RDTSC_START(comienzo);
                        Dt = getTraspuesta(D);
                        Dt_D = Dt * D;
                        RDTSC_STOP(final);
                        ciclos_clock[i] += final - comienzo;
                        ciclos_antes_de_leer_imagenes[i] = ciclos_clock[i];
                    }
                    for(size_t ind_arch = 0; ind_arch < archivos.size(); ++ind_arch){
                        vector<vector<double> > *imagen_entera = leerCSV(archivos[ind_arch]);
                        vector<vector<double> > imagen_discreta;
                        vector<double> vec_imagen_discreta;
                        vector<double> t_sin_ruido;
                        unsigned long ciclos_antes_del_ruido[repeticiones];
                        for(uint8_t i = 0; i < repeticiones; ++i) {
                            ciclos_clock[i] = ciclos_antes_de_leer_imagenes[i];
                            RDTSC_START(comienzo);
                            imagen_discreta = discretizar(*imagen_entera, discretizaciones[ind_disc]);
                            vec_imagen_discreta = pasarAVector(imagen_discreta);
                            t_sin_ruido = D * vec_imagen_discreta;
                            RDTSC_STOP(final);
                            ciclos_clock[i] += final - comienzo;
                            ciclos_antes_del_ruido[i] = ciclos_clock[i];
                        }
                        salida.open(nombre_arch_salida, ios::app);
                        salida << "Imagen: "+archivos[ind_arch] << endl;
                        salida.close();
                        for(size_t ind_ruido = 0; ind_ruido < ruidos.size(); ++ind_ruido){
                            double error;
                            salida.open(nombre_arch_salida, ios::app);
                            salida << "\t"<< "Ruido N°" << ind_ruido+1 << ":" << endl;
                            salida << "\t\tCiclos de clock: ";
                            for(uint8_t i = 0; i < repeticiones; ++i) {
                                ciclos_clock[i] = ciclos_antes_del_ruido[i];
                                RDTSC_START(comienzo);
                                vector<double> t_con_ruido = uniformNoise(t_sin_ruido, ruidos[ind_ruido].first, ruidos[ind_ruido].second, 0);
                                pair<vector<double>, short> v = EG2(Dt_D, Dt * t_con_ruido);
                                error = ECM(vec_imagen_discreta, v.first);
                                RDTSC_STOP(final);
                                ciclos_clock[i] += final - comienzo;
                                salida << "N°" << i+1 << " " << ciclos_clock[i] << "; ";
                            }
                            salida << endl;
                            salida << "\t\tError: " << error << endl;
                            salida.close();
                        }
                        salida.open(nombre_arch_salida, ios::app);
                        salida << endl << endl; //Lo hago 2 veces para mejor visibilidad.
                        salida.close();
                    }
                }
            }
        }
    }
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
    for (int j=0; j<ancho; j++) {
        for (int i = 0; i < ancho; i++) {
            valor = floor(vector[i + j]);
            linea += to_string((signed short) valor) + " ";
        }
        salida << linea << endl;
        linea = "";
    }
    salida.close();
}
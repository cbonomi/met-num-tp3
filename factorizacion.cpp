#include "factorizacion.h"

#define TOLERANCIA_REDONDEO 0.00001

using namespace std;


vector<double> sumaVec(const vector<double> &vec1, const vector<double> &vec2) {
    vector<double> res (vec1.size());
    for (uint i = 0; i < vec1.size(); i++)
        res[i] = vec1[i]+vec2[i];
    return res;
}

vector<double> restaVec(const vector<double> &vec1, const vector<double> &vec2) {
    vector<double> res (vec1.size());
    for (uint i = 0; i < vec1.size(); i++)
        res[i] = vec1[i]-vec2[i];
    return res;
}

double norma1(const vector<double> &v){
    double res = 0;
    for(size_t i = 0; i < v.size(); ++i)
        res += abs(v[i]);
    return res;
}

double norma2(const vector<double> &vec){//no tomo raiz porque no hace falta en nuestro caso
    double acum = 0;
    for (uint i = 0; i < vec.size(); i++) {
        acum+= vec[i]*vec[i];
    }
    return acum;
}

vector<double> mult_matr_por_vect(const vector<vector<double> > &M, const vector<double> &v){//recordar que el vector v son las inversas de las velocidades
    const size_t& cant_filas = M.size();
    const size_t& cant_columnas = M[0].size();
    if (cant_columnas != v.size()){
        cerr << "Los tamaños no coinciden " << cant_filas << "x" << cant_columnas << " * " << v.size() << "x1" <<  endl; 
    }
    vector<double> res(cant_filas);
    for(size_t i = 0; i < cant_filas; ++i) {
        res[i] = 0;
        for (size_t k = 0; k < cant_columnas; ++k)
            res[i] += M[i][k]*v[k];
    }
    return res;
}

void normalizar2(vector<double>& v){     //Según norma 2
    double norma = sqrt(norma2(v));
    for(size_t i = 0; i < v.size(); ++i)
        v[i] = v[i]/norma;
}

double producto_interno(const vector<double> &v, const vector<double> &v2){
    double res = 0;
    if(v.size() != v2.size())
        cout << "producto_interno: los vectores no son del mismo tamaño" << endl;
    else{
        for(size_t i = 0; i < v.size(); ++i)
            res += v[i]*v2[i];
    }
    return res;
}

pair<double,vector<double> > metodoPotencia(const vector<vector<double> > &M) {
    const size_t& n = M[0].size();
    double MAX_DIF = 0.0001;
    uint MAX_ITER = 500;
    pair<double,vector<double> > res2;
    double& autovalor2 = res2.first;
    double autovalor2_temp;
    vector<double>& autovector2 = res2.second;
    autovector2 = vector<double>(n,0);
    autovector2[0] = 1;   //Empieza siendo e1.
    vector<double> autovector2_temp;
    //Cálculo del autovalor:
    double diferencia2 = 1;
    float cantidad_iteraciones2 = 0;
    //auto t3 = chrono::system_clock::now();
    while(diferencia2 >= MAX_DIF and not (cantidad_iteraciones2 >= MAX_ITER)){
        autovector2_temp = mult_matr_por_vect(M, autovector2);
        autovalor2_temp = producto_interno(autovector2, autovector2_temp); //autovector está normalizado.
        normalizar2(autovector2_temp);
        autovector2 = mult_matr_por_vect(M, autovector2_temp);
        autovalor2 = producto_interno(autovector2_temp, autovector2); //autovector_temp está normalizado.
        normalizar2(autovector2);
        diferencia2 = abs(autovalor2-autovalor2_temp);
        ++cantidad_iteraciones2;
    }
    //auto t4 = chrono::system_clock::now();
    if(autovalor2 < 0){
        autovector2 = mult_matr_por_vect(M, autovector2);
        normalizar2(autovector2);
        diferencia2 = norma1(restaVec(autovector2, autovector2_temp));
        while(diferencia2 >= MAX_DIF and not (cantidad_iteraciones2 >= MAX_ITER)){
            autovector2_temp = mult_matr_por_vect(M, autovector2);
            normalizar2(autovector2_temp);
            autovector2 = mult_matr_por_vect(M, autovector2_temp);
            autovector2 = mult_matr_por_vect(M, autovector2);
            normalizar2(autovector2);
            diferencia2 = norma1(restaVec(autovector2, autovector2_temp));
            ++cantidad_iteraciones2;
        }
    }else{
        diferencia2 = norma1(restaVec(autovector2, autovector2_temp));
        while(diferencia2 >= MAX_DIF and not (cantidad_iteraciones2 >= MAX_ITER)){
            autovector2_temp = mult_matr_por_vect(M, autovector2);
            normalizar2(autovector2_temp);
            autovector2 = mult_matr_por_vect(M, autovector2_temp);
            normalizar2(autovector2);
            diferencia2 = norma1(restaVec(autovector2, autovector2_temp));
            ++cantidad_iteraciones2;
        }
    }
    return res2;
}


void multVecEsc(vector<double> &vec, double escalar) {//para multiplicar un vector por un escalar. Afecta al vector parámetro.
    for (uint i = 0; i<vec.size();i++)
            vec[i] *= escalar;
}
void multMatEsc(vector<vector<double> > &mat, double escalar) {//para multiplicar una matriz por un escalar. Afecta a la matriz parámetro.
    for (uint i = 0; i<mat.size();i++)
        for (uint j = 0; j<mat[i].size();j++)
            mat[i][j] *= escalar;
}

vector<vector<double> > multVec(const vector<double> &vec1) {//para generar una matriz a partir de un vector y su traspuesto
    const size_t& n = vec1.size();
    vector<vector<double> > res(n, vector<double>(n));
    for (uint i = 0; i<n;i++)
        for (uint j = 0; j<n;j++)
            res[i][j] = vec1[i]*vec1[j];
    return res;
}

vector<vector<double> > calcularXtX (const vector<vector<double> >& X) {
    const size_t& n = X.size();
    const size_t& m = X[0].size();
    vector<vector<double> > res(m, vector<double>(m));
    for (uint i = 0; i < m; ++i){
        for (uint j = i; j < m; ++j){
            for (uint k = 0; k < n; ++k){
                // Using X_T won't trash the cache.
                res[i][j] += X[k][i] * X[k][j];

            }
            res[j][i] = res[i][j];
        }
    }
    return res;
}

void sumMat(vector<vector<double> > &mat1, const vector<vector<double> > &mat2) {//suma de matrices
    for (uint i = 0; i<mat1.size();i++)
        for (uint j = 0; j<mat1.size();j++)
            mat1[i][j] += mat2[i][j];
}

vector< pair<double,vector<double> > > deflacion(vector<vector<double> > mat,const double tolerance, const bool debug) {
    uint N = mat[0].size();

    vector< pair<double,vector<double> > > res;    
    res.push_back(metodoPotencia(mat));
    vector<vector<double> > v_x_vt = multVec(res[0].second);    //v*vt
    multMatEsc(v_x_vt,res[0].first*(-1));                       //-lambda_i*(v*vt)
    sumMat(mat, v_x_vt);   
    for (uint i = 1; i < N; i++){  
        auto autov = metodoPotencia(mat);
        if (autov.first/res[0].first < tolerance){
                break;
        }
        res.push_back(autov);
        if (debug){
            cerr << "calculando autovalores: " << i+1 << "/" << N 
                 << " O_max:" << res[0].first << " O_i:" << res[i].first
                 << " Numero cond:" << res[0].first/res[i].first//'\r';
                 << " K+1/O1:" << res[i].first/res[0].first << endl;//'\r';
        }
        vector<vector<double> > v_x_vt = multVec(res[i].second);    //v*vt
        multMatEsc(v_x_vt,res[i].first*(-1));                       //-lambda_i*(v*vt)
        sumMat(mat, v_x_vt);
    }
    cout << "Numero de condicion: " << res[0].first/res[res.size()-1].first << endl;
    return res;
}

void calcular_svd(const vector<vector<double> > &mat,
                        vector<vector<double> > &Ut,    //
                        vector<double> &Sigm,
                        vector<vector<double> > &Vt,
                        const double tolerance,
                        const bool debug)
{
    uint N = mat.size();
    uint M = mat[0].size();
    if (debug)
        cerr << "Calculando SVD N=" << N << " M=" << M << endl;
    vector<vector<double> > AtA = calcularXtX(mat);
    vector< pair<double,vector<double> > > autovectores = deflacion(AtA,tolerance,debug);
    

    //Calculamos Sigma
    uint c = 0;
    Sigm.clear();
    while ((c<autovectores.size()) && (sqrt(autovectores[c].first) > TOLERANCIA_REDONDEO)){
        Sigm.push_back(sqrt(autovectores[c].first));
        c++;
    }

    if (debug)
        cerr << "TERMINO Sig tam=" << Sigm.size() << endl;
    //Calculamos Vt
    Vt.clear();
    for (uint j = 0; j < c; ++j){
            Vt.push_back(autovectores[j].second);
    }

    if (debug)
        cerr << "TERMINO Vt N=" << Vt.size() << " M=" << Vt[0].size() << endl;
    //Calculamos Ut
    Ut.clear();
    for (uint i = 0; i < Sigm.size(); ++i){
        vector<double> vec = Vt[i];
        multVecEsc(vec,double(1/Sigm[i]));
        vector<double> u_i = mult_matr_por_vect(mat,vec);
        
        Ut.push_back(u_i);
    }
    //Completo las filas de Ut
    if (debug)
        cerr << "TERMINO Ut N=" << Ut.size() << " M=" << Ut[0].size() << endl;
}
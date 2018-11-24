#ifndef PAGERANK_VECTORMAPMATRIX_H
#define PAGERANK_VECTORMAPMATRIX_H

#include <map>
#include <vector>

using namespace std;

class VectorMapMatrix {
public:
    VectorMapMatrix(); //Construyo una "matriz" de 0x0

    VectorMapMatrix(uint h, uint w); //Nueva matriz "Llena de ceros" de altura h, ancho w.

    VectorMapMatrix(const VectorMapMatrix &orig) = default; //default copy constructor

    VectorMapMatrix& operator=(const VectorMapMatrix &orig) = default; //default operador de asignacion

    ~VectorMapMatrix() = default; //destructor por defecto

    size_t cantFilas() const;

    size_t cantColumnas() const;

    void asignar(uint f, uint c, const double value); //Si el valor a asignar puede o no ser 0, usar esta función (y no operator[]).

    double at(uint f, uint c) const; //útil si queres leer la posición sin asignar un 0 (operator[] crea el nodo sin importar si no asignas nada).

    double& operator[](pair<uint,uint> p); //Usar solo si se quieren hacer muchas asignaciones distintas de 0.
    // Cuidado, no usar para asignar ceros, usar asignar en tal caso.

	map<uint, double>& operator[](uint i){ return m[i]; }; //Devuelve una referencia a la i-esima fila.

    VectorMapMatrix operator+(VectorMapMatrix const &B);

    vector<vector<double>> operator*(const VectorMapMatrix &B);

    vector<double> operator*(const vector<double> &v);

    void agregarFila(const map<uint,double>& map); //Agrega una nueva fila.

    void reservar(uint ancho, uint alto); //Permite reservar mas espacio para Maps si se desean agregar mas, y agrandar el ancho de los maps.
    //Solo lo cambia si es mayor, si es menor no cambia nada.

    //VectorMapMatrix mult(VectorMapMatrix const &B);

    //vector<double> mult(vector<double>&);

    void operator*(double valor);

    //VectorMapMatrix triangularMatriz();

	pair<vector<double>,short> EG(const VectorMapMatrix& thisTranspuestaAux, vector<double> b);    //Resolución de sist. de ec. mediante eliminación gaussiana (solo para matrices cuadradas).

	pair<vector<double>,short> EGPP(vector<double>);    //Resolución de sist. de ec. mediante eliminación gaussiana con pivoteo parcial (solo para matrices cuadradas).
	
	//VectorMapMatrix permutar(unsigned int j, unsigned int i);
private:
    vector<map<uint, double> > m; //La matriz va a tener un vector vertical, de arboles rojo negro horizontales.
    size_t width;

    void resta_de_filas(uint fila_a_modificar, double escalar, uint fila_para_restar);
};

ostream& operator << (ostream &o, VectorMapMatrix &B);

void mostrar_matriz_por_consola(VectorMapMatrix& m, string nombre_de_la_matriz);

void mostrar_vector_por_consola(vector<double>& v, string nombre_del_vector);

VectorMapMatrix vector2matrix(vector<double>& v, uint cant_filas);

#endif //PAGERANK_VECTORMAPMATRIX_H

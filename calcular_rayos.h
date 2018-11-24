#ifndef CALCULAR_RAYOS_H
#define CALCULAR_RAYOS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <functional>
#include <tuple>
#include <stdlib.h>
#include "VectorMapMatrix.h"
#include "util.h"


#include <cmath>

using namespace std;


/**
 * Dado 2 puntos por los que pasa una recta (que esten en el borde), los mueve para que la recta se rote en dirección
 * contraria a las agujas del reloj. Importante, estoy considerando que el pixel [0,0] será la esquina superior
 * izquierda, se debería cambiar esto si lo vamos a pensar distinto.
 * @param p1 es el punto 1 donde el primer valor es la fila, el segundo la columna, y debe estar en el borde de la
 * imagen
 * @param p2 es el segundo punto donde termina la recta, con las mismas condiciones que p1.
 * @param n es el tamaño total de la imágen final resultante de todos los métodos.
 */
void rotarContrarreloj(pair<uint,uint>& p1, pair<uint,uint>& p2, size_t n);

/**
 * Crea un par, de dos vectores (solo para tenerlos todos juntos y ordenados) donde el primer vector tiene en la
 * posición i-ésima al punto del pixel de inicio, de la recta rayo, que pasa por la fila i-ésima, el segundo vector
 * tiene el punto del pixel donde termina dicha recta.
 * @param n tamaño de la matriz imagen resultante, que será del mismo tamaño que la matriz D^k.
 * @param cada_cuanto distancia entre cada laser generado.
 * @param offset distancia inicial al primer laser.
 * @result par de vectores con los pixeles inicio fin de las rectas horizontales de toda una matriz imagen de nxn.
 */
pair<vector<pair<uint,uint> >, vector<pair<uint,uint> > > inicios_fines_horizontales(size_t n, size_t cada_cuanto, size_t offset);

/**
 * Crea una matriz D_k del mismo tamaño que la imágen que será el resultado (nxn), donde un 1 en esa posición indica
 * que la recta pasa por ese pixel. (osea, un 1 en D[0][37] indica que el rayo pasa por el pixel de la fila 1,
 * de la columna 38).
 * @param p1 es el primer punto, la función funciona para cualquier punto inicial EN EL BORDE
 * (no me fije que pasa si no esta en el borde, quizá ande)
 * A diferencia de matlab, decidi que para mejorar la comprensión, queda mas claro si lo pensamos matricial y no por ejes.
 * O sea, p1= [3,7], significa que p1 es el punto en el pixel de la fila 4, la columna 8. En matlab hice que sea
 * el punto 3 en el eje x, 7 en el eje y.
 * Por último el algoritmo supone que p1 y p2 son distintos, si son el mismo no funciona.
 * @param p2 es el segundo punto, donde termina la recta (lo mismo que el primer punto, puede ser cualquiera pero en borde)
 * @param n es el tamaño de la matriz resultado del método que estamos desarrollando (Por ahora sería cuadrada).
 * @return la matriz resultante con unos donde pasa el rayo, en este caso del tp sería la matriz D^k.
 */
vector<vector<double> > trazar_recta_en_matriz_D(pair<uint,uint> p1, pair<uint,uint> p2, size_t n);


/**
 * Carga las posiciones donde se colocarían laseres a lo alto (en el borde izquierdo y derecho) para una imagen de nxn.
 * esto indica la posicion de inicio de los laseres, su punto al que apuntan será manejado con otra función.
 * @param n: Tamaño de la matriz imágen donde se colocaran los laseres.
 * @param cada_cuanto: Cada cuantos pixeles se añadirá un nuevo laser, el primero se añade sin dejar espacio (si offset
 * es 0).
 * @param offset: espacio dejado antes del primer laser (se empiezan a colocar de arriba a abajo, solo en los bordes
 * izquierdo y derecho)
 * @param cant_maxima: parámetro extra para evitar que el método cree mas laseres que cant_maxima (por cada lado, o sea
 * se crean cant_maxima de laseres por el lado izquierdo, y cant_maxima por el lado derecho). Si cant_maxima es igual
 * a 0, se ignora el valor y se coloca tantos como sean posibles.
 * @return Devuelve los puntos donde inician los laseres para una matriz de nxn.
 */

vector<pair<uint,uint> > crearLaseres(size_t n, size_t cada_cuanto, size_t offset, size_t cant_maxima);


/**
 * Dado un vector con los laseres de inicio, crea los n puntos a donde empezarían a apuntar para dar una sola pasada.
 * @param Laseres: resultado dado por la función crearLaseres.
 * @param n: tamaño de la imágen.
 * @return vector con los puntos a los que apuntan los laseres, ordenados de la misma forma que se encuentran los
 * Laseres en el vector de Laseres.
 */
vector<pair<uint,uint> > crearPuntosDeFin(vector<pair<uint,uint> > Laseres, size_t n);


/**
 * Función que rota el punto donde termina el laser creado en crearPuntosDeFin. A los laseres que hayan empezado del
 * lado derecho, los rota contrarreloj (empezando de la esquina arriba a la derecha), y los laseres del lado izquierdo
 * los rota en dirección de las agujas del reloj, en el único caso en el que los laseres quedan horizontales, para no
 * repetir laseres, decidimos que los laseres de la izquierda roten un paso mas que los laseres de la derecha para que
 * no tengamos 2 veces laseres horizontales.
 * @param Laseres: vector de puntos creado con la función crearLaseres.
 * @param A_donde_apuntan: vector de puntos creado con la función crearPuntosDeFin, este lo toma por referencia y los
 * rota acorde a lo mencionado mas arriba.
 */

void barrerLaseres_H(const vector<pair<uint,uint> >& Laseres, vector<pair<uint,uint> >& A_donde_apuntan, size_t n);

/**
 * Igual que arriba pero no se saltean los horizontales, (esta incorporado en otra funcion el salto, ya que existen
 * otros repetidos, y queriamos evitarlos.
 */

void barrerLaseres_H_sin_salto(const vector<pair<uint,uint> >& Laseres, vector<pair<uint,uint> >& A_donde_apuntan, size_t n);

/**
 * Genera Matriz con todos los D_kij (cada fila es una de las matrices D_k).
 * @param tamMatriz tamaño de la imagen discretizada.
 * @param metodo_usado es un numero QUE DEBE VALER 0,1 o 2, y que indica, si es 0, que se usara el metodo de rotaciones
 * iniciando con rayos horizontales, si vale 1, serán unos rayos fijos, que son colocados en los lados horizontales de
 * la imagen y rotaran, si vale 2, estos rayos son colocados en el tope y fondo verticales de la imagen, y tambien rotan.
 * @param cantLaseres es la cantidad de laseres que se desean, DEBE SER DIVISOR DE tamMatriz o la función puede tener
 * resultados indeseables, (como minimo puede pasar que no se obtenga la cantidad deseada de laseres, o cosas peores).
 * @param saltear_hasta_n es la cantidad de pixeles rotados que saltearemos despues de cada rayo disparado, el minimo
 * valor permitido es 1 (CON 0 SE ROMPE) y aumentar el valor reduce el tiempo de computo, pero tambien reduce la
 * precision.
 * @return La matriz D con todos los D_k.
 */
VectorMapMatrix  generarRayos(size_t tamMatriz, int metodo_usado, int cantLaseres, int saltear_hasta_n);


/**
 * Esta funcion auxiliar sirve para chequear los rayos y puntos generados y poder graficar en matlab las rectas
 * correspondientes a los rayos generados por cada metodo, devuelve un par de 2 vectores, donde el primer vector guarda
 * los pares de enteros sin signo que representan los puntos de inicio de los rayos y el otro vector tendra los puntos
 * de finalizacion de dichos rayos para poder graficar en matlab.
 * @param tamMatriz tamaño de la imagen discretizada.
 * @param metodo_usado es un numero entero, y que indica, si es 0, que se usara el metodo de rotaciones
 * iniciando con rayos horizontales, si vale 1, serán unos rayos fijos, que son colocados en los lados horizontales de
 * la imagen y rotaran, si vale 2, estos rayos son colocados en el tope y fondo verticales de la imagen, y tambien rotan,
 * si vale 3 entonces se usa el metodo de horizontales agregando rayos que vengan del tope, si vale 4 usa un metodo en
 * el que evita repetir rayos, con rayos en la izquierda, derecha y arriba (aunque arriba hay algunos repetidos, 5 hace
 * sin repetidos solo horizontales (derecha e izquierda) y 6 hace sin repetidos verticales(rayos arriba y abajo) (para
 * cualquier valor distinto de los anteriores siempre vale el ultimo, aunque seria comportamiento no deseado).
 * @param cantLaseres es la cantidad de laseres que se desean, DEBE SER DIVISOR DE tamMatriz o la función puede tener
 * resultados indeseables, (como minimo puede pasar que no se obtenga la cantidad deseada de laseres, o cosas peores).
 * @param saltear_hasta_n es la cantidad de pixeles rotados que saltearemos despues de cada rayo disparado, el minimo
 * valor permitido es 1 (CON 0 SE ROMPE) y aumentar el valor reduce el tiempo de computo, pero tambien reduce la
 * precision.
 * @return La matriz D con todos los D_k.
 */
pair<vector<pair<uint,uint> >, vector<pair<uint,uint> > >  generarRayosDeControl(size_t tamMatriz, int metodo_usado, int cantLaseres, int saltear_hasta_n);

#endif //CALCULAR_RAYOS_H

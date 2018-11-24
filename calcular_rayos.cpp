#include "calcular_rayos.h"
/**
 * Dado 2 puntos por los que pasa una recta (que esten en el borde), los mueve para que la recta se rote en dirección
 * contraria a las agujas del reloj. Importante, estoy considerando que el pixel [0,0] será la esquina superior
 * izquierda, se debería cambiar esto si lo vamos a pensar distinto.
 * @param p1 es el punto 1 donde el primer valor es la fila, el segundo la columna, y debe estar en el borde de la
 * imagen
 * @param p2 es el segundo punto donde termina la recta, con las mismas condiciones que p1.
 * @param n es el tamaño total de la imágen final resultante de todos los métodos.
 */
void rotarContrarreloj(pair<uint,uint>& p1, pair<uint,uint>& p2, size_t n) {
    if(p1.second == 0 and p1.first != n-1) { //Si la columna es 0, entonces si la fila no es la última lo movemos para
        // abajo.
        p1.first++;
    } else if(p1.first == n-1 and p1.second != n-1) { //Estamos en última fila, y no estamos en última columna,
        // movemos a derecha.
        p1.second++;
    } else if(p1.second == n-1 and p1.first != 0) { //Última columna, y no en la primera fila, movemos para arriba.
        p1.first--;
    } else { //Estamos en la primera fila, y no estamos en la primera columna, entonces movemos a izquierda.
        p1.second--;
    }
    //repetimos lo anterior para p2.
    if(p2.second == 0 and p2.first != n-1) { //Si la columna es 0, entonces si la fila no es la última lo movemos para
        // abajo.
        p2.first++;
    } else if(p2.first == n-1 and p2.second != n-1) { //Estamos en última fila, y no estamos en última columna,
        // movemos a derecha.
        p2.second++;
    } else if(p2.second == n-1 and p2.first != 0) { //Última columna, y no en la primera fila, movemos para arriba.
        p2.first--;
    } else { //Estamos en la primera fila, y no estamos en la primera columna, entonces movemos a izquierda.
        p2.second--;
    }
    return;
}

/**
 * Crea un par, de dos vectores (solo para tenerlos todos juntos y ordenados) donde el primer vector tiene en la
 * posición i-ésima al punto del pixel de inicio, de la recta rayo, que pasa por la fila i-ésima, el segundo vector
 * tiene el punto del pixel donde termina dicha recta.
 * @param n tamaño de la matriz imagen resultante, que será del mismo tamaño que la matriz D^k.
 * @param cada_cuanto distancia entre cada laser generado, es deseable que cada_cuanto sea divisor de n, y offset
 * sea menor de cada_cuanto, si esto no se cumpliese puede que la funcion no haga toda la cantidad de laseres que se
 * deseaba.
 * @param offset distancia inicial al primer laser.
 * @result par de vectores con los pixeles inicio fin de las rectas horizontales de toda una matriz imagen de nxn.
 */
pair<vector<pair<uint,uint> >, vector<pair<uint,uint> > > inicios_fines_horizontales(size_t n, size_t cada_cuanto, size_t offset) {
    pair<vector<pair<uint,uint> >, vector<pair<uint,uint> > > result =
            make_pair(vector<pair<uint,uint> >(), vector<pair<uint,uint> >());
    uint columna_ini = 0;
    uint columna_fin = n-1;
    for(uint i = offset; i < n; i += cada_cuanto) {
        uint fila_ini = i;
        uint fila_fin = i;
        result.first.emplace_back(make_pair(fila_ini,columna_ini));
        result.second.emplace_back(make_pair(fila_fin,columna_fin));
    }
    return result;
}


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
vector<vector<double> > trazar_recta_en_matriz_D(pair<uint,uint> p1, pair<uint,uint> p2, size_t n) {
    vector<vector<double> > result(n, vector<double>(n,0));
    int inicio;
    int fin;
    if(p1.second == p2.second) {
        // La recta es vertical, la trazo y termino el método.
        for(uint i = 0; i<n ; i++) {
            result[i][p1.second] = 1;
        }
        return result;
    } else if(p1.second < p2.second) { //me fijo cual va primero.
        inicio = p1.second;
        fin = p2.second;
    } else {
        fin = p1.second;
        inicio = p2.second;
    }
    double a = (double(p2.first) - double(p1.first)) / (double(p2.second) - double(p1.second)); // No importa el orden
    // de p1 y p2, calcula la pendiente, fila 2 menos fila 1 sobre columna 2 menos columna 1.
    double b = double(p1.first) - a * double(p1.second) + 0.5; // calcula b = y - ax, con el punto p1.
    for(int i = inicio; i<fin ; i++) { //voy de la columna inicio a la fin, pintando los pixeles por los que pase.
        int pintar_desde = floor(a*i + b);
        int pintar_hasta = floor(a*(i+1) +b);
        if (pintar_hasta < pintar_desde) {
            int aux = pintar_desde; //si estan al revez los swappeo
            pintar_desde = pintar_hasta;
            pintar_hasta = aux;
        }
        for(int j = pintar_desde; j <= pintar_hasta ; j++) {
            result[j][i] = 1;
        }
    }
    int pintar_desde = floor(a*fin + b);
    int pintar_hasta = floor(a*(fin+1) +b);
    if(pintar_hasta >= int(n)) {
        pintar_hasta = n-1;
    } else if (pintar_hasta < 0) {
        pintar_hasta = 0;
    }
    if (pintar_hasta < pintar_desde) {
        int aux = pintar_desde; //si estan al revez los swappeo
        pintar_desde = pintar_hasta;
        pintar_hasta = aux;
    }

    for(int j = pintar_desde; j <= pintar_hasta ; j++) {
        result[j][fin] = 1;
    }
    return result;
}

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

vector<pair<uint,uint> > crearLaseres(size_t n, size_t cada_cuanto, size_t offset, size_t cant_maxima){
    vector<pair<uint,uint> > result;
    uint cant_creada = 0;
    for(uint i = offset; i<n and (cant_maxima == 0 or cant_creada < cant_maxima); i += cada_cuanto) {
        result.emplace_back(make_pair(i,uint(0)));
	    ++cant_creada;
    }

    cant_creada = 0;

    for(uint i = offset; i<n and (cant_maxima == 0 or cant_creada < cant_maxima); i += cada_cuanto) {
        result.emplace_back(i,uint(n-1));
	    ++cant_creada;
    }

    return result;
}


/**
 * Dado un vector con los laseres de inicio, crea los n puntos a donde empezarían a apuntar para dar una sola pasada.
 * @param Laseres: resultado dado por la función crearLaseres.
 * @param n: tamaño de la imágen.
 * @return vector con los puntos a los que apuntan los laseres, ordenados de la misma forma que se encuentran los
 * Laseres en el vector de Laseres.
 */

vector<pair<uint,uint> > crearPuntosDeFin(vector<pair<uint,uint> > Laseres, size_t n) {
    vector<pair<uint,uint> > result;
    result.reserve(Laseres.size()*2);
    for(uint i = 0; i<Laseres.size(); i++) {
        uint fila;
        uint columna;
        if(Laseres[i].second == 0){
            fila = uint(0);
            columna = uint(1);
        } else {
            fila = uint(0);
            columna = uint(n-2);
        }
        pair<uint,uint> p = make_pair(fila,columna);
        result.emplace_back(p);
    }
    return result;
}


/**
 * Función que rota el punto donde termina el laser creado en crearPuntosDeFin. A los laseres que hayan empezado del
 * lado derecho, los rota contrarreloj (empezando de la esquina arriba a la derecha), y los laseres del lado izquierdo
 * los rota en dirección de las agujas del reloj, en el unico caso en el que los laseres quedan horizontales, para no
 * repetir laseres, decidimos que los laseres de la izquierda roten un paso más, y asi saltearse el valor horizontal,
 * ya que los laseres de la derecha ya lo representan, para que no tengamos 2 veces laseres horizontales.
 * @param Laseres: vector de puntos creado con la función crearLaseres.
 * @param A_donde_apuntan: vector de puntos creado con la función crearPuntosDeFin, este lo toma por referencia y los
 * rota acorde a lo mencionado mas arriba.
 */

void barrerLaseres_H(const vector<pair<uint,uint> >& Laseres, vector<pair<uint,uint> >& A_donde_apuntan, size_t n) {
    for(uint i = 0; i<Laseres.size(); i++) {
        if(Laseres[i].second == 0){ //roto en direccion del reloj.
            if(A_donde_apuntan[i].second == 0 and A_donde_apuntan[i].first != 0) { //Si la columna es 0, entonces si la fila no es la primera lo movemos para
                // arriba.
                A_donde_apuntan[i].first--;
            } else if(A_donde_apuntan[i].first == n-1 and A_donde_apuntan[i].second != 0) { //Estamos en última fila, y no estamos en primera columna,
                // movemos a izquierda.
                A_donde_apuntan[i].second--;
            } else if(A_donde_apuntan[i].second == n-1 and A_donde_apuntan[i].first != n-1) { //Última columna, y no en la ultima fila, movemos para abajo.
                A_donde_apuntan[i].first++;
                if (A_donde_apuntan[i].first == Laseres[i].first and A_donde_apuntan[i].first != n-1) { //Si el laser de la izquierda queda horizontal, lo muevo hacia abajo, excepto si el laser es el ultimo mas abajo.
                    A_donde_apuntan[i].first++;
                } else if (A_donde_apuntan[i].first == Laseres[i].first and A_donde_apuntan[i].first == n-1) { // Si es el ultimo mas abajo, lo muevo hacia la izquierda (porque sino se va de rango).
                    A_donde_apuntan[i].second--;
                }
            } else { //Estamos en la primera fila, y no estamos en la última columna, entonces movemos a derecha.
                A_donde_apuntan[i].second++;
                if(Laseres[i].first == 0 and A_donde_apuntan[i].second == n-1) { // Si es el primer laser y cae en la ultima columna, entonces lo debemos mover hacia abajo ya que no lo queremos horizontal.
                    A_donde_apuntan[i].first++;
                }
            }
        } else { //rota contrarreloj.
            if(A_donde_apuntan[i].second == 0 and A_donde_apuntan[i].first != n-1) { //Si la columna es 0, entonces si la fila no es la última lo movemos para
                // abajo.
                A_donde_apuntan[i].first++;
            } else if(A_donde_apuntan[i].first == n-1 and A_donde_apuntan[i].second != n-1) { //Estamos en última fila, y no estamos en última columna,
                // movemos a derecha.
                A_donde_apuntan[i].second++;
            } else if(A_donde_apuntan[i].second == n-1 and A_donde_apuntan[i].first != 0) { //Última columna, y no en la primera fila, movemos para arriba.
                A_donde_apuntan[i].first--;
            } else { //Estamos en la primera fila, y no estamos en la primera columna, entonces movemos a izquierda.
                A_donde_apuntan[i].second--;
            }
        }
    }
    return;
}

void barrerLaseres_H_sin_salto(const vector<pair<uint,uint> >& Laseres, vector<pair<uint,uint> >& A_donde_apuntan, size_t n) {
    for(uint i = 0; i<Laseres.size(); i++) {
        if(Laseres[i].second == 0){ //roto en direccion del reloj.
            if(A_donde_apuntan[i].second == 0 and A_donde_apuntan[i].first != 0) { //Si la columna es 0, entonces si la fila no es la primera lo movemos para
                // arriba.
                A_donde_apuntan[i].first--;
            } else if(A_donde_apuntan[i].first == n-1 and A_donde_apuntan[i].second != 0) { //Estamos en última fila, y no estamos en primera columna,
                // movemos a izquierda.
                A_donde_apuntan[i].second--;
            } else if(A_donde_apuntan[i].second == n-1 and A_donde_apuntan[i].first != n-1) { //Última columna, y no en la ultima fila, movemos para abajo.
                A_donde_apuntan[i].first++;
            } else { //Estamos en la primera fila, y no estamos en la última columna, entonces movemos a derecha.
                A_donde_apuntan[i].second++;
            }
        } else { //rota contrarreloj.
            if(A_donde_apuntan[i].second == 0 and A_donde_apuntan[i].first != n-1) { //Si la columna es 0, entonces si la fila no es la última lo movemos para
                // abajo.
                A_donde_apuntan[i].first++;
            } else if(A_donde_apuntan[i].first == n-1 and A_donde_apuntan[i].second != n-1) { //Estamos en última fila, y no estamos en última columna,
                // movemos a derecha.
                A_donde_apuntan[i].second++;
            } else if(A_donde_apuntan[i].second == n-1 and A_donde_apuntan[i].first != 0) { //Última columna, y no en la primera fila, movemos para arriba.
                A_donde_apuntan[i].first--;
            } else { //Estamos en la primera fila, y no estamos en la primera columna, entonces movemos a izquierda.
                A_donde_apuntan[i].second--;
            }
        }
    }
    return;
}

/**
 * Genera Matriz con todos los D_kij (cada fila es una de las matrices D_k).
 * @param tamMatriz tamaño de la imagen discretizada.
 * @param metodo_usado es un numero entero que indica:
 *      Si vale 0, que se usara el metodo de rotaciones iniciando con rayos horizontales.
 *      Si vale 1, serán unos rayos fijos, que son colocados en los lados horizontales de la imagen y rotaran.
 *      Si vale 2, estos rayos son colocados en el tope y fondo verticales de la imagen, y tambien rotan.
 *      Si vale 3, entonces se usa el metodo de horizontales agregando rayos que vengan del tope.
 *      Si vale 4, usa un metodo en el que evita repetir rayos, con rayos en la izquierda, derecha y arriba
 *                 (aunque arriba hay algunos repetidos).
 *      Si vale 5, hace sin repetidos solo horizontales (derecha e izquierda)
 *      Si vale 6, hace sin repetidos verticales (rayos arriba y abajo)
 * (para cualquier valor distinto de los anteriores siempre vale el ultimo, aunque seria comportamiento no deseado).
 * @param cantLaseres es la cantidad de laseres que se desean, DEBE SER DIVISOR DE tamMatriz o la función puede tener
 * resultados indeseables, (como minimo puede pasar que no se obtenga la cantidad deseada de laseres, o cosas peores).
 * @param saltear_hasta_n es la cantidad de pixeles rotados que saltearemos despues de cada rayo disparado, el minimo
 * valor permitido es 1 (CON 0 SE ROMPE) y aumentar el valor reduce el tiempo de computo, pero tambien reduce la
 * precision.
 * @return La matriz D con todos los D_k.
 */
VectorMapMatrix  generarRayos(size_t tamMatriz, int metodo_usado, int cantLaseres, int saltear_hasta_n) {
    // creamos un laser de cada una de las esquinas
    if (metodo_usado == 0){
        pair<vector<pair<uint,uint> >, vector<pair<uint,uint> > > laseresYsensores =
                inicios_fines_horizontales(tamMatriz, tamMatriz/cantLaseres, tamMatriz/(cantLaseres*2));

        vector<pair<uint,uint> > laseres = laseresYsensores.first;
        vector<pair<uint,uint> > sensores = laseresYsensores.second;

        VectorMapMatrix D_ks(0, tamMatriz*tamMatriz); /*creamos el result con 0 maps (los vamos a agregar despues
 * uno por uno.*/
        D_ks.reservar(tamMatriz*tamMatriz, 2 * cantLaseres * tamMatriz); /* Este es el vector con las matrices D, para cada uno de los K rayos (hay que convertirlas en vectores).
        Tenemos cantLaseres rayos que rotaremos aproximadamente 2n veces (y asi los rotamos 180º). */

        vector<vector<double> > D_k; //matriz auxiliar del D_k del laser a calcular.
        int rotaciones = 0; //Cantidad de rotaciones ejecutadas, solo calcularemos cuando rotaciones%saltear_hasta_n == 0
        for(uint i = 0; i < 2 * tamMatriz; i++) {
            if (rotaciones % saltear_hasta_n == saltear_hasta_n/2){ //si rotamos la cantidad correcta entonces calculamos.
                for(uint j = 0; j < laseres.size(); j++) {
                    D_k = trazar_recta_en_matriz_D(laseres[j], sensores[j], tamMatriz);
                    map<uint, double> D_k_map = pasarAMap(D_k);
                    D_ks.agregarFila(D_k_map);
                }
            }

            for(uint j = 0; j<laseres.size(); j++) {
                rotarContrarreloj(laseres[j],sensores[j],tamMatriz);
            }
            rotaciones++;
        }

        return D_ks;


    } else if (metodo_usado == 1) {
        vector<pair<uint,uint> > laseres = crearLaseres(tamMatriz, tamMatriz/cantLaseres, tamMatriz/(cantLaseres*2), 0); //tamano, despues cada_cuanta_dist, offset, max_cant de rayos.
        vector<pair<uint,uint> > sensores = crearPuntosDeFin(laseres, tamMatriz);

        VectorMapMatrix D_ks(0, tamMatriz*tamMatriz);

        D_ks.reservar(tamMatriz*tamMatriz, 6 * tamMatriz*cantLaseres);
        /* Este es el vector con las matrices D, para cada uno de los K rayos (hay que convertirlas en vectores).
        Tenemos 2*cantLaseres rayos que rotaremos aproximadamente 3tamMatriz veces. */

        vector<vector<double> > D_k; //matriz auxiliar del D_k del laser a calcular.

        int rotaciones = 0; //Cantidad de rotaciones ejecutadas, solo calcularemos cuando rotaciones%saltear_hasta_n == 0
        while(sensores[0].first != tamMatriz - 1 or sensores[0].second != 0) { //Esto quiza es dificil de ver, pero para los laseres izquierdos
            // que se saltean el horizontal, este es la ultima posicion interesante a la que apuntan. NOTA IMPORTANTE,
            // SI SE HACEN MAS DE UN SALTO PUEDE QUE ESTO NO TERMINE. ASIQUE CUIDADO CON PONER MAS DE UN rotarLaseres.
            if (rotaciones % saltear_hasta_n == saltear_hasta_n/2) { //si rotamos la cantidad correcta entonces calculamos.
                for (uint i = 0; i < laseres.size(); i++) {
                    D_k = trazar_recta_en_matriz_D(laseres[i], sensores[i], tamMatriz);
                    map<uint, double> D_k_map = pasarAMap(D_k);
                    D_ks.agregarFila(D_k_map);
                }
            }
            barrerLaseres_H(laseres,sensores,tamMatriz);
            rotaciones++;
        }

        return D_ks;


    } else if (metodo_usado == 2) {
        vector<pair<uint,uint> > laseres = crearLaseres(tamMatriz, tamMatriz/cantLaseres, tamMatriz/(cantLaseres*2), 0); //tamano, despues cada_cuanta_dist, offset, max_cant de rayos.
        vector<pair<uint,uint> > sensores = crearPuntosDeFin(laseres, tamMatriz);

        VectorMapMatrix D_ks(0, tamMatriz*tamMatriz);

        D_ks.reservar(tamMatriz*tamMatriz, 6 * tamMatriz*cantLaseres); //ancho, alto
        /* Este es el vector con las matrices D, para cada uno de los K rayos (hay que convertirlas en vectores).
        Tenemos 2*cantLaseres rayos que rotaremos aproximadamente 3tamMatriz veces. */

        vector<vector<double> > D_k; //matriz auxiliar del D_k del laser a calcular.
        int rotaciones = 0; //Cantidad de rotaciones ejecutadas, solo calcularemos cuando rotaciones%saltear_hasta_n == 0


        while(sensores[0].first != tamMatriz - 1 or sensores[0].second != 0) { //Esto quiza es dificil de ver, pero para los laseres izquierdos
            // que se saltean el horizontal, este es la ultima posicion interesante a la que apuntan. NOTA IMPORTANTE,
            // SI SE HACEN MAS DE UN SALTO PUEDE QUE ESTO NO TERMINE. ASIQUE CUIDADO CON PONER MAS DE UN rotarLaseres.
            if (rotaciones % saltear_hasta_n == saltear_hasta_n /2) { //si rotamos la cantidad correcta entonces calculamos.
                for(uint i = 0; i < laseres.size(); i++) {
                    D_k = trazar_recta_en_matriz_D(laseres[i], sensores[i], tamMatriz);
                    //Traspongo para conseguir la recta si el rayo fuese vertical.
                    vector<vector<double> > D_k_transp(D_k.size(), vector<double>(D_k[0].size()) );

                    for(uint i = 0; i < D_k.size(); ++i)
                        for (unsigned int j=0; j < D_k[0].size(); ++j)
                            D_k_transp[j][i] = D_k[i][j];


                    map<uint, double> D_k_map = pasarAMap(D_k_transp);
                    D_ks.agregarFila(D_k_map);
                }
            }
            barrerLaseres_H(laseres,sensores,tamMatriz);
            rotaciones++;
        }

        return D_ks;
    } else if (metodo_usado == 3) {
        vector<pair<uint,uint> > laseres = crearLaseres(tamMatriz, tamMatriz/cantLaseres, tamMatriz/(cantLaseres*2), 0); //tamano, despues cada_cuanta_dist, offset, max_cant de rayos.
        vector<pair<uint,uint> > sensores = crearPuntosDeFin(laseres, tamMatriz);

        VectorMapMatrix D_ks(0, tamMatriz*tamMatriz);

        D_ks.reservar(tamMatriz*tamMatriz, 6 * tamMatriz*cantLaseres); //ancho, alto
        /* Este es el vector con las matrices D, para cada uno de los K rayos (hay que convertirlas en vectores).
        Tenemos 2*cantLaseres rayos que rotaremos aproximadamente 3tamMatriz veces. */

        vector<vector<double> > D_k; //matriz auxiliar del D_k del laser a calcular.
        int rotaciones = 0; //Cantidad de rotaciones ejecutadas, solo calcularemos cuando rotaciones%saltear_hasta_n == 0


        while(sensores[0].first != tamMatriz - 1 or sensores[0].second != 0) { //Esto quiza es dificil de ver, pero para los laseres izquierdos
            // que se saltean el horizontal, este es la ultima posicion interesante a la que apuntan. NOTA IMPORTANTE,
            // SI SE HACEN MAS DE UN SALTO PUEDE QUE ESTO NO TERMINE. ASIQUE CUIDADO CON PONER MAS DE UN rotarLaseres.
            if (rotaciones % saltear_hasta_n == saltear_hasta_n/2) { //si rotamos la cantidad correcta entonces calculamos.
                for(uint i = 0; i < laseres.size(); i++) {
                    D_k = trazar_recta_en_matriz_D(laseres[i], sensores[i], tamMatriz);
                    //Traspongo para conseguir la recta si el rayo fuese vertical.
                    vector<vector<double> > D_k_transp(D_k.size(), vector<double>(D_k[0].size()) );

                    map<uint, double> D_k_map = pasarAMap(D_k);
                    D_ks.agregarFila(D_k_map);

                    if (i>=laseres.size()/2) {
                        if(sensores[i].second == 0) { //calculo el transpuesto que es cambiar fila y columna para los verticales.
                            D_k = trazar_recta_en_matriz_D(make_pair(laseres[i].second,laseres[i].first),
                                                           make_pair(sensores[i].second,sensores[i].first), tamMatriz);
                            map<uint, double> D_k_map = pasarAMap(D_k);
                            D_ks.agregarFila(D_k_map);
                        }
                    }
                }
            }
            barrerLaseres_H(laseres,sensores,tamMatriz);
            rotaciones++;
        }

        return D_ks;
    } else if(metodo_usado == 4) {
        vector<pair<uint,uint> > laseres = crearLaseres(tamMatriz, tamMatriz/cantLaseres, tamMatriz/(cantLaseres*2), 0); //tamano, despues cada_cuanta_dist, offset, max_cant de rayos.
        vector<pair<uint,uint> > sensores = crearPuntosDeFin(laseres, tamMatriz);

        vector<pair<uint,uint> >::iterator it_ini_laseres = laseres.begin();
        vector<pair<uint,uint> >::iterator it_med_laseres = laseres.begin() + laseres.size()/2;
        vector<pair<uint,uint> >::iterator it_fin_laseres = laseres.end();

        vector<pair<uint,uint> >::iterator it_ini_sensores = sensores.begin();
        vector<pair<uint,uint> >::iterator it_med_sensores = sensores.begin() + sensores.size()/2;
        vector<pair<uint,uint> >::iterator it_fin_sensores = sensores.end();

        vector<pair<uint,uint> > laseres0(it_ini_laseres, it_med_laseres); //laseres izq
        vector<pair<uint,uint> > laseres1(it_med_laseres, it_fin_laseres); //laseres der
        vector<pair<uint,uint> > sensores0(it_ini_sensores, it_med_sensores); //sensores izq
        vector<pair<uint,uint> > sensores1(it_med_sensores, it_fin_sensores); //sensores der

        VectorMapMatrix D_ks(0, tamMatriz*tamMatriz);

        D_ks.reservar(tamMatriz*tamMatriz, 6 * tamMatriz*cantLaseres); //ancho, alto
        /* Este es el vector con las matrices D, para cada uno de los K rayos (hay que convertirlas en vectores).
        Tenemos 2*cantLaseres rayos que rotaremos aproximadamente 3tamMatriz veces. */

        vector<vector<double> > D_k; //matriz auxiliar del D_k del laser a calcular.
        int rotaciones = 0; //Cantidad de rotaciones ejecutadas, solo calcularemos cuando rotaciones%saltear_hasta_n == 0


        while(sensores1[0].first != tamMatriz - 1 or sensores1[0].second != tamMatriz - 1) { //Esto quiza es dificil de ver, pero para los laseres derechos
            // esta es la ultima posicion interesante a la que apuntan. NOTA IMPORTANTE,
            // SI SE HACEN MAS DE UN SALTO PUEDE QUE ESTO NO TERMINE. ASIQUE CUIDADO CON PONER MAS DE UN rotarLaseres.
            if (rotaciones % saltear_hasta_n == saltear_hasta_n /2) { //si rotamos la cantidad correcta entonces calculamos.
                for(uint i = 0; i < laseres1.size(); i++) {
                    D_k = trazar_recta_en_matriz_D(laseres1[i], sensores1[i], tamMatriz);
                    map<uint, double> D_k_map = pasarAMap(D_k);
                    D_ks.agregarFila(D_k_map);
                    if(sensores1[i].second == 0) { //calculo el transpuesto que es cambiar fila y columna para los verticales.
                        D_k = trazar_recta_en_matriz_D(make_pair(laseres1[i].second,laseres1[i].first),
                                make_pair(sensores1[i].second,sensores1[i].first), tamMatriz);
                        map<uint, double> D_k_map = pasarAMap(D_k);
                        D_ks.agregarFila(D_k_map);
                    }
                }
            }
            barrerLaseres_H_sin_salto(laseres1,sensores1,tamMatriz);
            rotaciones++;
        }

        rotaciones = 0;

        while(sensores0[0].first != tamMatriz - 1 or sensores0[0].second != 0) { //Esto quiza es dificil de ver, pero para los laseres izquierdos
            // este es la ultima posicion interesante a la que apuntan. NOTA IMPORTANTE,
            // SI SE HACEN MAS DE UN SALTO PUEDE QUE ESTO NO TERMINE. ASIQUE CUIDADO CON PONER MAS DE UN rotarLaseres.
            if (rotaciones % saltear_hasta_n == saltear_hasta_n /2) { //si rotamos la cantidad correcta entonces calculamos.
                bool chequearSensores = true; /*Quiero chequear si hay un laser en laseres1 repetido con el sensor actual
 * para evitar repetidos */
                for (uint i = 0; i < laseres1.size(); i++) {
                    if(laseres1[i] == sensores0[0]) { //todos los sensores apuntan al mismo lugar, por lo que elegimos el primero
                        chequearSensores = false; // Encontramos coincidencia, entonces no calculamos D_k.
                    }
                }
                if (chequearSensores) { //Si chequeamos que no habia laser en la posicion, calculamos las D_k.
                    for(uint i = 0; i < laseres0.size(); i++) {
                        D_k = trazar_recta_en_matriz_D(laseres0[i], sensores0[i], tamMatriz);
                        map<uint, double> D_k_map = pasarAMap(D_k);
                        D_ks.agregarFila(D_k_map);
                    }
                }
            }
            barrerLaseres_H_sin_salto(laseres0,sensores0,tamMatriz);
            rotaciones++;
        }

        return D_ks;
    } else if (metodo_usado == 5){
        vector<pair<uint,uint> > laseres = crearLaseres(tamMatriz, tamMatriz/cantLaseres, tamMatriz/(cantLaseres*2), 0); //tamano, despues cada_cuanta_dist, offset, max_cant de rayos.
        vector<pair<uint,uint> > sensores = crearPuntosDeFin(laseres, tamMatriz);

        vector<pair<uint,uint> >::iterator it_ini_laseres = laseres.begin();
        vector<pair<uint,uint> >::iterator it_med_laseres = laseres.begin() + laseres.size()/2;
        vector<pair<uint,uint> >::iterator it_fin_laseres = laseres.end();

        vector<pair<uint,uint> >::iterator it_ini_sensores = sensores.begin();
        vector<pair<uint,uint> >::iterator it_med_sensores = sensores.begin() + sensores.size()/2;
        vector<pair<uint,uint> >::iterator it_fin_sensores = sensores.end();

        vector<pair<uint,uint> > laseres0(it_ini_laseres, it_med_laseres); //laseres izq
        vector<pair<uint,uint> > laseres1(it_med_laseres, it_fin_laseres); //laseres der
        vector<pair<uint,uint> > sensores0(it_ini_sensores, it_med_sensores); //sensores izq
        vector<pair<uint,uint> > sensores1(it_med_sensores, it_fin_sensores); //sensores der

        VectorMapMatrix D_ks(0, tamMatriz*tamMatriz);

        D_ks.reservar(tamMatriz*tamMatriz, 6 * tamMatriz*cantLaseres); //ancho, alto
        /* Este es el vector con las matrices D, para cada uno de los K rayos (hay que convertirlas en vectores).
        Tenemos 2*cantLaseres rayos que rotaremos aproximadamente 3tamMatriz veces. */

        vector<vector<double> > D_k; //matriz auxiliar del D_k del laser a calcular.
        int rotaciones = 0; //Cantidad de rotaciones ejecutadas, solo calcularemos cuando rotaciones%saltear_hasta_n == 0


        while(sensores1[0].first != tamMatriz - 1 or sensores1[0].second != tamMatriz - 1) { //Esto quiza es dificil de ver, pero para los laseres derechos
            // esta es la ultima posicion interesante a la que apuntan. NOTA IMPORTANTE,
            // SI SE HACEN MAS DE UN SALTO PUEDE QUE ESTO NO TERMINE. ASIQUE CUIDADO CON PONER MAS DE UN rotarLaseres.
            if (rotaciones % saltear_hasta_n == saltear_hasta_n /2) { //si rotamos la cantidad correcta entonces calculamos.
                for(uint i = 0; i < laseres1.size(); i++) {
                    D_k = trazar_recta_en_matriz_D(laseres1[i], sensores1[i], tamMatriz);
                    map<uint, double> D_k_map = pasarAMap(D_k);
                    D_ks.agregarFila(D_k_map);
                }
            }
            barrerLaseres_H_sin_salto(laseres1,sensores1,tamMatriz);
            rotaciones++;
        }

        rotaciones = 0;

        while(sensores0[0].first != tamMatriz - 1 or sensores0[0].second != 0) { //Esto quiza es dificil de ver, pero para los laseres izquierdos
            // este es la ultima posicion interesante a la que apuntan. NOTA IMPORTANTE,
            // SI SE HACEN MAS DE UN SALTO PUEDE QUE ESTO NO TERMINE. ASIQUE CUIDADO CON PONER MAS DE UN rotarLaseres.
            if (rotaciones % saltear_hasta_n == saltear_hasta_n /2) { //si rotamos la cantidad correcta entonces calculamos.
                bool chequearSensores = true; /*Quiero chequear si hay un laser en laseres1 repetido con el sensor actual
 * para evitar repetidos */
                for (uint i = 0; i < laseres1.size(); i++) {
                    if(laseres1[i] == sensores0[0]) { //todos los sensores apuntan al mismo lugar, por lo que elegimos el primero
                        chequearSensores = false; // Encontramos coincidencia, entonces no calculamos D_k.
                    }
                }
                if (chequearSensores) { //Si chequeamos que no habia laser en la posicion, calculamos las D_k.
                    for(uint i = 0; i < laseres0.size(); i++) {
                        D_k = trazar_recta_en_matriz_D(laseres0[i], sensores0[i], tamMatriz);
                        map<uint, double> D_k_map = pasarAMap(D_k);
                        D_ks.agregarFila(D_k_map);
                    }
                }
            }
            barrerLaseres_H_sin_salto(laseres0,sensores0,tamMatriz);
            rotaciones++;
        }

        return D_ks;
    } else {
        vector<pair<uint,uint> > laseres = crearLaseres(tamMatriz, tamMatriz/cantLaseres, tamMatriz/(cantLaseres*2), 0); //tamano, despues cada_cuanta_dist, offset, max_cant de rayos.
        vector<pair<uint,uint> > sensores = crearPuntosDeFin(laseres, tamMatriz);

        vector<pair<uint,uint> >::iterator it_ini_laseres = laseres.begin();
        vector<pair<uint,uint> >::iterator it_med_laseres = laseres.begin() + laseres.size()/2;
        vector<pair<uint,uint> >::iterator it_fin_laseres = laseres.end();

        vector<pair<uint,uint> >::iterator it_ini_sensores = sensores.begin();
        vector<pair<uint,uint> >::iterator it_med_sensores = sensores.begin() + sensores.size()/2;
        vector<pair<uint,uint> >::iterator it_fin_sensores = sensores.end();

        vector<pair<uint,uint> > laseres0(it_ini_laseres, it_med_laseres); //laseres izq
        vector<pair<uint,uint> > laseres1(it_med_laseres, it_fin_laseres); //laseres der
        vector<pair<uint,uint> > sensores0(it_ini_sensores, it_med_sensores); //sensores izq
        vector<pair<uint,uint> > sensores1(it_med_sensores, it_fin_sensores); //sensores der

        VectorMapMatrix D_ks(0, tamMatriz*tamMatriz);

        D_ks.reservar(tamMatriz*tamMatriz, 6 * tamMatriz*cantLaseres); //ancho, alto
        /* Este es el vector con las matrices D, para cada uno de los K rayos (hay que convertirlas en vectores).
        Tenemos 2*cantLaseres rayos que rotaremos aproximadamente 3tamMatriz veces. */

        vector<vector<double> > D_k; //matriz auxiliar del D_k del laser a calcular.
        int rotaciones = 0; //Cantidad de rotaciones ejecutadas, solo calcularemos cuando rotaciones%saltear_hasta_n == 0


        while(sensores1[0].first != tamMatriz - 1 or sensores1[0].second != tamMatriz - 1) { //Esto quiza es dificil de ver, pero para los laseres derechos
            // esta es la ultima posicion interesante a la que apuntan. NOTA IMPORTANTE,
            // SI SE HACEN MAS DE UN SALTO PUEDE QUE ESTO NO TERMINE. ASIQUE CUIDADO CON PONER MAS DE UN rotarLaseres.
            if (rotaciones % saltear_hasta_n == saltear_hasta_n /2) { //si rotamos la cantidad correcta entonces calculamos.
                for(uint i = 0; i < laseres1.size(); i++) {
                    D_k = trazar_recta_en_matriz_D(laseres1[i], sensores1[i], tamMatriz);
                    //Traspongo para conseguir la recta si el rayo fuese vertical.
                    vector<vector<double> > D_k_transp(D_k.size(), vector<double>(D_k[0].size()) );

                    for(uint i = 0; i < D_k.size(); ++i)
                        for (unsigned int j=0; j < D_k[0].size(); ++j)
                            D_k_transp[j][i] = D_k[i][j];


                    map<uint, double> D_k_map_trans = pasarAMap(D_k_transp);
                    D_ks.agregarFila(D_k_map_trans);
                }
            }
            barrerLaseres_H_sin_salto(laseres1,sensores1,tamMatriz);
            rotaciones++;
        }

        rotaciones = 0;

        while(sensores0[0].first != tamMatriz - 1 or sensores0[0].second != 0) { //Esto quiza es dificil de ver, pero para los laseres izquierdos
            // este es la ultima posicion interesante a la que apuntan. NOTA IMPORTANTE,
            // SI SE HACEN MAS DE UN SALTO PUEDE QUE ESTO NO TERMINE. ASIQUE CUIDADO CON PONER MAS DE UN rotarLaseres.
            if (rotaciones % saltear_hasta_n == saltear_hasta_n /2) { //si rotamos la cantidad correcta entonces calculamos.
                bool chequearSensores = true; /*Quiero chequear si hay un laser en laseres1 repetido con el sensor actual
 * para evitar repetidos */
                for (uint i = 0; i < laseres1.size(); i++) {
                    if(laseres1[i] == sensores0[0]) { //todos los sensores apuntan al mismo lugar, por lo que elegimos el primero
                        chequearSensores = false; // Encontramos coincidencia, entonces no calculamos D_k.
                    }
                }
                if (chequearSensores) { //Si chequeamos que no habia laser en la posicion, calculamos las D_k.
                    for(uint i = 0; i < laseres0.size(); i++) {
                        D_k = trazar_recta_en_matriz_D(laseres0[i], sensores0[i], tamMatriz);
                        //Traspongo para conseguir la recta si el rayo fuese vertical.
                        vector<vector<double> > D_k_transp(D_k.size(), vector<double>(D_k[0].size()) );

                        for(uint i = 0; i < D_k.size(); ++i)
                            for (unsigned int j=0; j < D_k[0].size(); ++j)
                                D_k_transp[j][i] = D_k[i][j];


                        map<uint, double> D_k_map_trans = pasarAMap(D_k_transp);
                        D_ks.agregarFila(D_k_map_trans);
                    }
                }
            }
            barrerLaseres_H_sin_salto(laseres0,sensores0,tamMatriz);
            rotaciones++;
        }

        return D_ks;
    }
}

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
pair<vector<pair<uint,uint> >, vector<pair<uint,uint> > >  generarRayosDeControl(size_t tamMatriz, int metodo_usado, int cantLaseres, int saltear_hasta_n) {
    // creamos un laser de cada una de las esquinas
    if (metodo_usado == 0){
        pair<vector<pair<uint,uint> >, vector<pair<uint,uint> > > laseresYsensores =
                inicios_fines_horizontales(tamMatriz, tamMatriz/cantLaseres, tamMatriz/(cantLaseres*2));

        vector<pair<uint,uint> > laseres = laseresYsensores.first;
        vector<pair<uint,uint> > sensores = laseresYsensores.second;

        pair<vector<pair<uint,uint> >, vector<pair<uint,uint> > > result;

        int rotaciones = 0; //Cantidad de rotaciones ejecutadas, solo calcularemos cuando rotaciones%saltear_hasta_n == 0
        for(uint i = 0; i < 2 * tamMatriz; i++) {
            if (rotaciones % saltear_hasta_n == saltear_hasta_n/2){ //si rotamos la cantidad correcta entonces calculamos.
                for(uint j = 0; j < laseres.size(); j++) {
                    result.first.emplace_back(laseres[j]);
                    result.second.emplace_back(sensores[j]);
                }
            }

            for(uint j = 0; j<laseres.size(); j++) {
                rotarContrarreloj(laseres[j],sensores[j],tamMatriz);
            }
            rotaciones++;
        }

        return result;


    } else if (metodo_usado == 1) {
        vector<pair<uint,uint> > laseres = crearLaseres(tamMatriz, tamMatriz/cantLaseres, tamMatriz/(cantLaseres*2), 0); //tamano, despues cada_cuanta_dist, offset, max_cant de rayos.
        vector<pair<uint,uint> > sensores = crearPuntosDeFin(laseres, tamMatriz);

        pair<vector<pair<uint,uint> >, vector<pair<uint,uint> > > result;

        int rotaciones = 0; //Cantidad de rotaciones ejecutadas, solo calcularemos cuando rotaciones%saltear_hasta_n == 0
        while(sensores[0].first != tamMatriz - 1 or sensores[0].second != 0) { //Esto quiza es dificil de ver, pero para los laseres izquierdos
            // que se saltean el horizontal, este es la ultima posicion interesante a la que apuntan. NOTA IMPORTANTE,
            // SI SE HACEN MAS DE UN SALTO PUEDE QUE ESTO NO TERMINE. ASIQUE CUIDADO CON PONER MAS DE UN rotarLaseres.
            if (rotaciones % saltear_hasta_n == saltear_hasta_n/2) { //si rotamos la cantidad correcta entonces calculamos.
                for (uint i = 0; i < laseres.size(); i++) {
                    result.first.emplace_back(laseres[i]);
                    result.second.emplace_back(sensores[i]);
                }
            }
            barrerLaseres_H(laseres,sensores,tamMatriz);
            rotaciones++;
        }

        return result;


    } else if (metodo_usado == 2) {
        vector<pair<uint,uint> > laseres = crearLaseres(tamMatriz, tamMatriz/cantLaseres, tamMatriz/(cantLaseres*2), 0); //tamano, despues cada_cuanta_dist, offset, max_cant de rayos.
        vector<pair<uint,uint> > sensores = crearPuntosDeFin(laseres, tamMatriz);

        pair<vector<pair<uint,uint> >, vector<pair<uint,uint> > > result;

        int rotaciones = 0; //Cantidad de rotaciones ejecutadas, solo calcularemos cuando rotaciones%saltear_hasta_n == 0


        while(sensores[0].first != tamMatriz - 1 or sensores[0].second != 0) { //Esto quiza es dificil de ver, pero para los laseres izquierdos
            // que se saltean el horizontal, este es la ultima posicion interesante a la que apuntan. NOTA IMPORTANTE,
            // SI SE HACEN MAS DE UN SALTO PUEDE QUE ESTO NO TERMINE. ASIQUE CUIDADO CON PONER MAS DE UN rotarLaseres.
            if (rotaciones % saltear_hasta_n == saltear_hasta_n /2) { //si rotamos la cantidad correcta entonces calculamos.
                for(uint i = 0; i < laseres.size(); i++) {
                    result.first.emplace_back(make_pair(laseres[i].second, laseres[i].first));
                    result.second.emplace_back(make_pair(sensores[i].second,sensores[i].first));
                }
            }
            barrerLaseres_H(laseres,sensores,tamMatriz);
            rotaciones++;
        }

        return result;

    } else if (metodo_usado == 3) {
        vector<pair<uint,uint> > laseres = crearLaseres(tamMatriz, tamMatriz/cantLaseres, tamMatriz/(cantLaseres*2), 0); //tamano, despues cada_cuanta_dist, offset, max_cant de rayos.
        vector<pair<uint,uint> > sensores = crearPuntosDeFin(laseres, tamMatriz);

        pair<vector<pair<uint,uint> >, vector<pair<uint,uint> > > result;

        int rotaciones = 0; //Cantidad de rotaciones ejecutadas, solo calcularemos cuando rotaciones%saltear_hasta_n == 0


        while(sensores[0].first != tamMatriz - 1 or sensores[0].second != 0) { //Esto quiza es dificil de ver, pero para los laseres izquierdos
            // que se saltean el horizontal, este es la ultima posicion interesante a la que apuntan. NOTA IMPORTANTE,
            // SI SE HACEN MAS DE UN SALTO PUEDE QUE ESTO NO TERMINE. ASIQUE CUIDADO CON PONER MAS DE UN rotarLaseres.
            if (rotaciones % saltear_hasta_n == saltear_hasta_n/2) { //si rotamos la cantidad correcta entonces calculamos.
                for(uint i = 0; i < laseres.size(); i++) {
                    result.first.emplace_back(laseres[i]);
                    result.second.emplace_back(sensores[i]);

                    if (i>=laseres.size()/2) {
                        if(sensores[i].second == 0) { //si son rayos horizontales de derecha los copio verticales.
                            result.first.emplace_back(make_pair(laseres[i].second, laseres[i].first));
                            result.second.emplace_back(make_pair(sensores[i].second,sensores[i].first));
                        }
                    }
                }
            }
            barrerLaseres_H(laseres,sensores,tamMatriz);
            rotaciones++;
        }

        return result;
    } else if(metodo_usado == 4) {
        vector<pair<uint,uint> > laseres = crearLaseres(tamMatriz, tamMatriz/cantLaseres, tamMatriz/(cantLaseres*2), 0); //tamano, despues cada_cuanta_dist, offset, max_cant de rayos.
        vector<pair<uint,uint> > sensores = crearPuntosDeFin(laseres, tamMatriz);

        vector<pair<uint,uint> >::iterator it_ini_laseres = laseres.begin();
        vector<pair<uint,uint> >::iterator it_med_laseres = laseres.begin() + laseres.size()/2;
        vector<pair<uint,uint> >::iterator it_fin_laseres = laseres.end();

        vector<pair<uint,uint> >::iterator it_ini_sensores = sensores.begin();
        vector<pair<uint,uint> >::iterator it_med_sensores = sensores.begin() + sensores.size()/2;
        vector<pair<uint,uint> >::iterator it_fin_sensores = sensores.end();

        vector<pair<uint,uint> > laseres0(it_ini_laseres, it_med_laseres); //laseres izq
        vector<pair<uint,uint> > laseres1(it_med_laseres, it_fin_laseres); //laseres der
        vector<pair<uint,uint> > sensores0(it_ini_sensores, it_med_sensores); //sensores izq
        vector<pair<uint,uint> > sensores1(it_med_sensores, it_fin_sensores); //sensores der

        pair<vector<pair<uint,uint> >, vector<pair<uint,uint> > > result;

        int rotaciones = 0; //Cantidad de rotaciones ejecutadas, solo calcularemos cuando rotaciones%saltear_hasta_n == 0


        while(sensores1[0].first != tamMatriz - 1 or sensores1[0].second != tamMatriz - 1) { //Esto quiza es dificil de ver, pero para los laseres derechos
            // esta es la ultima posicion interesante a la que apuntan. NOTA IMPORTANTE,
            // SI SE HACEN MAS DE UN SALTO PUEDE QUE ESTO NO TERMINE. ASIQUE CUIDADO CON PONER MAS DE UN rotarLaseres.
            if (rotaciones % saltear_hasta_n == saltear_hasta_n /2) { //si rotamos la cantidad correcta entonces calculamos.
                for(uint i = 0; i < laseres1.size(); i++) {
                    result.first.emplace_back(laseres1[i]);
                    result.second.emplace_back(sensores1[i]);

                    if(sensores1[i].second == 0) { //si estan dando al lado opuesto, tambien coloco los transpuestos.
                        result.first.emplace_back(make_pair(laseres1[i].second, laseres1[i].first));
                        result.second.emplace_back(make_pair(sensores1[i].second,sensores1[i].first));
                    }
                }
            }
            barrerLaseres_H_sin_salto(laseres1,sensores1,tamMatriz);
            rotaciones++;
        }

        rotaciones = 0;

        while(sensores0[0].first != tamMatriz - 1 or sensores0[0].second != 0) { //Esto quiza es dificil de ver, pero para los laseres izquierdos
            // este es la ultima posicion interesante a la que apuntan. NOTA IMPORTANTE,
            // SI SE HACEN MAS DE UN SALTO PUEDE QUE ESTO NO TERMINE. ASIQUE CUIDADO CON PONER MAS DE UN rotarLaseres.
            if (rotaciones % saltear_hasta_n == saltear_hasta_n /2) { //si rotamos la cantidad correcta entonces calculamos.
                bool chequearSensores = true; /*Quiero chequear si hay un laser en laseres1 repetido con el sensor actual
 * para evitar repetidos */
                for (uint i = 0; i < laseres1.size(); i++) {
                    if(laseres1[i] == sensores0[0]) { //todos los sensores apuntan al mismo lugar, por lo que elegimos el primero
                        chequearSensores = false; // Encontramos coincidencia, entonces no calculamos D_k.
                    }
                }
                if (chequearSensores) { //Si chequeamos que no habia laser en la posicion, calculamos las D_k.
                    for(uint i = 0; i < laseres0.size(); i++) {
                        result.first.emplace_back(laseres0[i]);
                        result.second.emplace_back(sensores0[i]);
                    }
                }
            }
            barrerLaseres_H_sin_salto(laseres0,sensores0,tamMatriz);
            rotaciones++;
        }

        return result;
    } else if (metodo_usado == 5){
        vector<pair<uint,uint> > laseres = crearLaseres(tamMatriz, tamMatriz/cantLaseres, tamMatriz/(cantLaseres*2), 0); //tamano, despues cada_cuanta_dist, offset, max_cant de rayos.
        vector<pair<uint,uint> > sensores = crearPuntosDeFin(laseres, tamMatriz);

        vector<pair<uint,uint> >::iterator it_ini_laseres = laseres.begin();
        vector<pair<uint,uint> >::iterator it_med_laseres = laseres.begin() + laseres.size()/2;
        vector<pair<uint,uint> >::iterator it_fin_laseres = laseres.end();

        vector<pair<uint,uint> >::iterator it_ini_sensores = sensores.begin();
        vector<pair<uint,uint> >::iterator it_med_sensores = sensores.begin() + sensores.size()/2;
        vector<pair<uint,uint> >::iterator it_fin_sensores = sensores.end();

        vector<pair<uint,uint> > laseres0(it_ini_laseres, it_med_laseres); //laseres izq
        vector<pair<uint,uint> > laseres1(it_med_laseres, it_fin_laseres); //laseres der
        vector<pair<uint,uint> > sensores0(it_ini_sensores, it_med_sensores); //sensores izq
        vector<pair<uint,uint> > sensores1(it_med_sensores, it_fin_sensores); //sensores der

        pair<vector<pair<uint,uint> >, vector<pair<uint,uint> > > result;

        int rotaciones = 0; //Cantidad de rotaciones ejecutadas, solo calcularemos cuando rotaciones%saltear_hasta_n == 0


        while(sensores1[0].first != tamMatriz - 1 or sensores1[0].second != tamMatriz - 1) { //Esto quiza es dificil de ver, pero para los laseres derechos
            // esta es la ultima posicion interesante a la que apuntan. NOTA IMPORTANTE,
            // SI SE HACEN MAS DE UN SALTO PUEDE QUE ESTO NO TERMINE. ASIQUE CUIDADO CON PONER MAS DE UN rotarLaseres.
            if (rotaciones % saltear_hasta_n == saltear_hasta_n /2) { //si rotamos la cantidad correcta entonces calculamos.
                for(uint i = 0; i < laseres1.size(); i++) {
                    result.first.emplace_back(laseres1[i]);
                    result.second.emplace_back(sensores1[i]);
                }
            }
            barrerLaseres_H_sin_salto(laseres1,sensores1,tamMatriz);
            rotaciones++;
        }

        rotaciones = 0;

        while(sensores0[0].first != tamMatriz - 1 or sensores0[0].second != 0) { //Esto quiza es dificil de ver, pero para los laseres izquierdos
            // este es la ultima posicion interesante a la que apuntan. NOTA IMPORTANTE,
            // SI SE HACEN MAS DE UN SALTO PUEDE QUE ESTO NO TERMINE. ASIQUE CUIDADO CON PONER MAS DE UN rotarLaseres.
            if (rotaciones % saltear_hasta_n == saltear_hasta_n /2) { //si rotamos la cantidad correcta entonces calculamos.
                bool chequearSensores = true; /*Quiero chequear si hay un laser en laseres1 repetido con el sensor actual
 * para evitar repetidos */
                for (uint i = 0; i < laseres1.size(); i++) {
                    if(laseres1[i] == sensores0[0]) { //todos los sensores apuntan al mismo lugar, por lo que elegimos el primero
                        chequearSensores = false; // Encontramos coincidencia, entonces no calculamos D_k.
                    }
                }
                if (chequearSensores) { //Si chequeamos que no habia laser en la posicion, calculamos las D_k.
                    for(uint i = 0; i < laseres0.size(); i++) {
                        result.first.emplace_back(laseres0[i]);
                        result.second.emplace_back(sensores0[i]);
                    }
                }
            }
            barrerLaseres_H_sin_salto(laseres0,sensores0,tamMatriz);
            rotaciones++;
        }

        return result;
    } else {
        vector<pair<uint,uint> > laseres = crearLaseres(tamMatriz, tamMatriz/cantLaseres, tamMatriz/(cantLaseres*2), 0); //tamano, despues cada_cuanta_dist, offset, max_cant de rayos.
        vector<pair<uint,uint> > sensores = crearPuntosDeFin(laseres, tamMatriz);

        vector<pair<uint,uint> >::iterator it_ini_laseres = laseres.begin();
        vector<pair<uint,uint> >::iterator it_med_laseres = laseres.begin() + laseres.size()/2;
        vector<pair<uint,uint> >::iterator it_fin_laseres = laseres.end();

        vector<pair<uint,uint> >::iterator it_ini_sensores = sensores.begin();
        vector<pair<uint,uint> >::iterator it_med_sensores = sensores.begin() + sensores.size()/2;
        vector<pair<uint,uint> >::iterator it_fin_sensores = sensores.end();

        vector<pair<uint,uint> > laseres0(it_ini_laseres, it_med_laseres); //laseres izq
        vector<pair<uint,uint> > laseres1(it_med_laseres, it_fin_laseres); //laseres der
        vector<pair<uint,uint> > sensores0(it_ini_sensores, it_med_sensores); //sensores izq
        vector<pair<uint,uint> > sensores1(it_med_sensores, it_fin_sensores); //sensores der

        pair<vector<pair<uint,uint> >, vector<pair<uint,uint> > > result;

        int rotaciones = 0; //Cantidad de rotaciones ejecutadas, solo calcularemos cuando rotaciones%saltear_hasta_n == 0


        while(sensores1[0].first != tamMatriz - 1 or sensores1[0].second != tamMatriz - 1) { //Esto quiza es dificil de ver, pero para los laseres derechos
            // esta es la ultima posicion interesante a la que apuntan. NOTA IMPORTANTE,
            // SI SE HACEN MAS DE UN SALTO PUEDE QUE ESTO NO TERMINE. ASIQUE CUIDADO CON PONER MAS DE UN rotarLaseres.
            if (rotaciones % saltear_hasta_n == saltear_hasta_n /2) { //si rotamos la cantidad correcta entonces calculamos.
                for(uint i = 0; i < laseres1.size(); i++) {
                    result.first.emplace_back(make_pair(laseres1[i].second, laseres1[i].first));
                    result.second.emplace_back(make_pair(sensores1[i].second,sensores1[i].first));
                }
            }
            barrerLaseres_H_sin_salto(laseres1,sensores1,tamMatriz);
            rotaciones++;
        }

        rotaciones = 0;

        while(sensores0[0].first != tamMatriz - 1 or sensores0[0].second != 0) { //Esto quiza es dificil de ver, pero para los laseres izquierdos
            // este es la ultima posicion interesante a la que apuntan. NOTA IMPORTANTE,
            // SI SE HACEN MAS DE UN SALTO PUEDE QUE ESTO NO TERMINE. ASIQUE CUIDADO CON PONER MAS DE UN rotarLaseres.
            if (rotaciones % saltear_hasta_n == saltear_hasta_n /2) { //si rotamos la cantidad correcta entonces calculamos.
                bool chequearSensores = true; /*Quiero chequear si hay un laser en laseres1 repetido con el sensor actual
 * para evitar repetidos */
                for (uint i = 0; i < laseres1.size(); i++) {
                    if(laseres1[i] == sensores0[0]) { //todos los sensores apuntan al mismo lugar, por lo que elegimos el primero
                        chequearSensores = false; // Encontramos coincidencia, entonces no calculamos D_k.
                    }
                }
                if (chequearSensores) { //Si chequeamos que no habia laser en la posicion, calculamos las D_k.
                    for(uint i = 0; i < laseres0.size(); i++) {
                        result.first.emplace_back(make_pair(laseres0[i].second, laseres0[i].first));
                        result.second.emplace_back(make_pair(sensores0[i].second,sensores0[i].first));
                    }
                }
            }
            barrerLaseres_H_sin_salto(laseres0,sensores0,tamMatriz);
            rotaciones++;
        }

        return result;
    }
}

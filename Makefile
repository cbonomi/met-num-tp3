CC=g++
CFLAGS= -O3 -std=c++11 -Wextra -Wall -pedantic -g -lm -Wno-unused-variable -Wno-unused-parameter



main: main.cpp util.h util.cpp calcular_rayos.h calcular_rayos.cpp VectorMapMatrix.h VectorMapMatrix.cpp factorizacion.h factorizacion.cpp
	$(CC) $(CFLAGS) $^ -o tp3


clean:
	rm tp2
	rm *.txt


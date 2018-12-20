
Instrucciones:
1. Convertir las imagenes a reconstruir con el script provisto por la catedra:

- Para convertir todas las imagenes con extensión "png" de una carpeta "imagenes/" en csv y poner las imagenes convertidas en el directorio "imagenes_convertidas/" correr:
    python csv_converter.py imagenes/ imagenes_convertidas/ .png


2. Compilar el codigo del TP de C++ con make.
3. Realizar la reconstruccion de la imagen con ./tp3 -i <ruta_input_csv> -o <ruta_output>. Tanto -i como -o son obligatorios para ejecutar el codigo. Si alguno no fue ingresado, se dara un mensaje de error indicando esto y se mostraran todos los posibles argumentos que utilizamos.


Ejemplo: ./tp3 -r 0.01 -t A -i tomo3.csv -o out.csv -nt 0.00 -m 1 -d 16
Genera la reconstruccion de tamaño 16 a partir de los datos en el archivo
  'tomo3.csv' utilizando el metodo de Barrido Horizontal con ruido aditivo 0.01


4. Para visualizar el resultado usar el script provisto por la catedra: csv_visualizer.py
Ejemplo: python csv_visualizer.py output.csv




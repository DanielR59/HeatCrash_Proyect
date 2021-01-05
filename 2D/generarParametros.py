import sys
from funciones import hdf5
if __name__ == "__main__":
    try:
        in_file_name=sys.argv[1]
    except:
        mensaje = """ Error: La ejecucion de este programa requiere de 1 argumento.
        Ejecucion correcta: python hdf5.py ENTRADA
        donde el nombre "ENTRADA" es el nombre del archivo donde
        se almacenan los datos del problema.
        Los datos que se guardan en el archivo ENTRADA son los que se encuentran en
        este archivo hdf5.py

        Por ejemplo: python hdf5.py ENTRADA"""

        print(mensaje)
        sys.exit(1)

    #Generamos parametros
    
    Datos={

        'ax' : 0, #Punto inicial en dirección x
        'bx' : 1,  #Punto final en dirección x
        'ay' : 0,  #Punto inicial en dirección y
        'by' : 1, #Punto final en dirección y
        'Nx' : 60, # Número de nodos en x
        'Ny' : 60, # Número de nodos en y
        'Tx1' : 40, # Condición de frontera 1 en x
        'Tx2' : -20, #Condición de frontera 2 en x
        'Ty1' : 3000, # Condición de frontera 1 en y
        'Ty2' : -100, # Condición de frontera 2 en x
        'Tini' : 120, # Condición inicial
        'ht' : 0.01, # Espaciamiento en tiempo
        'Tmax' : 3, # Tiempo máximo
        'Tolerancia' : 1E-4, # Tolerancia
        'kappa_x' : 1,
        'kappa_y' : 2,
        'fuente' : 0,
        'c_p' : 1,
        'rho' :2,
        'vel_x' : 1,
        'vel_y' : 1,
    }

    hdf5.saveParametros(in_file_name,Datos)

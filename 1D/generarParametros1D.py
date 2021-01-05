import sys
from Funciones1D import hdf5
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

        'a' : 0,
        'b': 2.5, 
        'N' : 300,
        'Ta': 1, 
        'Tb' : 0, 
        'k': 0.01, 
        'q': 0,
        'ht' : 0.001,
        'Tmax' : 1,
        'vel' : 1,
        'Tol' : 1E-4,
        'rho' : 1,
        'Nt' : 500,
        

    }

    hdf5.saveParametros(in_file_name,Datos)

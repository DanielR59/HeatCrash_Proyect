import hdf5
import numpy as np
import sys
import matplotlib.pyplot as plt


if __name__ == "__main__":
    
    try:
        in_file_name = sys.argv[1]; 
    except:
        mensaje = """ Error: La ejecucion de este programa requiere de 2 argumentos.
        Ejecucion correcta: python 2D_Poisson_01.py entrada salida
        donde "entrada" es el nombre de un archivo que contiene los
        datos del problema :  este se puede generar con el programa hdf5.py.
        El nombre "salida" se usa para almacenar la solucion_animada del problema.

        Por ejemplo: python 2D_Poisson_01.py SALIDA""" 

        print(mensaje)
        sys.exit(1)

    
    #leemos los parametros
    Datos=hdf5.leerParametros(in_file_name,'xg','yg','solucion','error')
    
    for key,val in Datos.items():

        exec(key + '=val')


    fig = plt.figure(figsize=(11, 7), dpi=100)
    plt.contourf(xg, yg, solucion)
    plt.colorbar()
    plt.contour(xg, yg, solucion)
    plt.xlabel('X')
    plt.ylabel('Y')
    
    try: 
        error
        plt.figure()
        plt.plot(error)
        plt.semilogy()
        
    except NameError:
        pass

    plt.show()

from funciones import hdf5
import numpy as np
import sys
import matplotlib.pyplot as plt


if __name__ == "__main__":
    
    try:
        in_file_name = sys.argv[1]; 
    except:
        mensaje = """ Error: La ejecucion de este programa requiere de 1 argumento.
        Ejecucion correcta: python {} ArchivoEntrada 

        Donde "ArchivoEntrada" es el nombre de un archivo que contiene la solucion almacenada 
        de un problema 2D

        Por ejemplo: python {} SALIDA""".format(__file__,__file__) 

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

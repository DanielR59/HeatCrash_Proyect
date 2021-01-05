from funciones import hdf5
from funciones.funciones2D import Animacion_Contorno, Animacion_Superficie
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation

if __name__ == "__main__":

    try:
        in_file_name = sys.argv[1];
    except:
        mensaje = """ Error: La ejecucion de este programa requiere de 2 argumentos.
        Ejecucion correcta: python {} entrada salida
        donde "entrada" es el nombre de un archivo que contiene los
        datos resultados de un problema 2D que involucre el tiempo.
        El nombre "salida" es el nombre y tipo de archivo que se generara de la animación.
        Este puede ser tanto .gif como .mp4
        NOTA: Si no tienes ffmpeg no podras guardar como .mp4

        Por ejemplo: python {} ENTRADA """.format(__file__,__file__)

        print(mensaje)
        sys.exit(1)

    #leemos los parametros
    Datos=hdf5.leerParametros(in_file_name,'Nt','ht','xg','yg','Zcambio')
    for key,val in Datos.items():
        exec(key + '=val')
    #Seleccionamos un paso adecuado
    Animacion_Contorno(xg,yg,Zcambio,Nt)
    Animacion_Superficie(xg,yg,Zcambio,Nt) 
    
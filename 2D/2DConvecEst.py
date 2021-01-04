from funciones.funciones2D import iterationConv2D, boundary_cond_dirichtlet, GrafContornos, Graf3D
from funciones import hdf5
import numpy as np
import sys
import time
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix

# =============================================================================
# Especificación de apertura correcta del programa
# =============================================================================

if __name__ == "__main__":
    #                       Tx2
    #  0,0                 frontera A
    #   o-----o-----o-----o-----o-----o-----o-----o
    #   |     |     |     |     |     |     |     |     f
    #   o-----o-----o-----o-----o-----o-----o-----o     r
    #   |     |     |     |     |     |     |     |     o
    # C o-----o-----o-----o-----o-----o-----o-----o     n
    # C |     |     |     |     |     |     |     | Ny  t
    # C o-----o-----o-----o-----o-----o-----o-----o     e
    #   |     |     |     |     |       |   |     |
    #   o-----o-----o-----o-----o-----o-----o-----o    D
    #   |     |     |     |     |     |     |     |
    #   o-----o-----o-----o-----o-----o-----o-----o
    #                       Nx
    #                       Tx1

    #                   frontera B
    try:
        in_file_name = sys.argv[1]; out_file_name = sys.argv[2]
    except:
        mensaje = """ Error: La ejecucion de este programa requiere de 2 argumentos.
        Ejecucion correcta: python {} entrada salida
        donde "entrada" es el nombre de un archivo que contiene los
        datos del problema :  este se puede generar con el programa hdf5.py.
        El nombre "salida" se usa para almacenar la solucion del problema.

        Por ejemplo: python {} ENTRADA SALIDA""".format(__file__,__file__)

        print(mensaje)
        sys.exit(1)

# =============================================================================
# Lectura de parámetros del archivo de entrada
# =============================================================================
    Datos=hdf5.leerParametros(in_file_name,'ax','ay','bx','by','kappa','Nx','Ny','Tx1','Tx2','Ty1','Ty2','Tolerancia','fuente','c_p','rho','vel_x','vel_y','kappa_x','kappa_y')
    
# =============================================================================
# Imprimimos los parametros y los evaluamos
# =============================================================================
    for key,val in Datos.items(): #ciclo evaluador
        print(key,'=',val)
        print('-'*10)
        exec(key + '=val')

    hx = (bx-ax)/(Nx+1) # espaciamiento entre nodos dirección X
    hy = (by-ay)/(Ny+1) # espaciamiento entre nodos dirección Y
    x = np.linspace(ax,bx,Nx+2) # dominio en dirección X
    y = np.linspace(ay,by,Ny+2) # dominio en dirección Y
    xg, yg = np.meshgrid(x,y) # generación de malla


    print('hx = ',hx)
    print('hy = ',hy)
# =============================================================================
# Definición de la mátriz del sistema
# =============================================================================
    u = np.zeros((Ny+2, Nx+2))
    #Aplicación de condiciones de frontera tipo dirichlet
    u= boundary_cond_dirichtlet(u,Tx1,Tx2,Ty1,Ty2)
    #Este es el lado derecho de la ecuación, que contiene la condición inicial
    f = np.ones_like(u)*fuente # RHS

    # Definicion de constantes para el calculo
    alfa_x=c_p*rho*vel_x
    alfa_y=c_p*rho*vel_y
# =============================================================================
# Definición las condiciones de Tolerancia: criterio de término anticipado
# =============================================================================    

    for i in range(20000):
        u,error=iterationConv2D(u,f,alfa_x,alfa_y,kappa_x,kappa_y,hx,hy)
        if error < Tolerancia:
            break
# =============================================================================
# Conservar los datos resultantes del calculo en formato hdf5 
# =============================================================================      
    Datos['xg']=xg
    Datos['yg']=yg
    Datos['solucion']=u

    hdf5.saveParametros(out_file_name,Datos)
# =============================================================================
# Generación de la gráfica
# =============================================================================
    #figura 1: mapa de contornos
    GrafContornos(xg, yg, u, 8, 0.75, 'inferno')
    
    #figura 2: superficie
    Graf3D(xg, yg, u, 'inferno')
        

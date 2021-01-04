# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 13:15:40 2020

@author: Kelly
"""

import numpy as np
import matplotlib.pyplot as plt
from Funciones1D import funcionesCE1D as func
from Funciones1D import hdf5 
import sys
if __name__ == "__main__":
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


    #Se leen los parámetros de interés del archivo generado por el programa hdf51D.py desde terminal
    par =hdf5.leerParametros(in_file_name,'a','b','rho', 'N', 'k', 'vel', 'Ta', 'Tb' ) 
    print('Los parámetros empleados son:',par)

    a = par['a']
    b = par['b']
    L = b-a
    N = par['N']
    h= L/(N+1)
    print('El valor de h es:', h)
    k = par['k']
    rho = par['rho']
    vel = par['vel']
    Ta = par['Ta']
    Tb = par['Tb']

    #Se manda a llamar la función para contruir a matriz A
    A=func.Laplaciano1DConvEst(N, h, k, rho, vel)

    #Se imprime la matriz A
    print('La matriz A es:',A)

    # Aplicacion de las cond. Dirichlet para el calculo correcto de la sol.
    f= np.zeros(N)
    f[0] = Ta*((rho*vel)/(2*h) + (k)/(h**2)) 
    f[N-1]= -Tb*((rho*vel)/(2*h) - (k)/(h**2))
    print('',f)


    # La solucion en el arreglo u, incluye las fronteras
    U = np.zeros(N+2)
    #obtención de solucion del sistema de N x N

    U[1:N+1] = np.linalg.solve(A,f)

    #Valores de frontera conocidos debido a las condiciones de Dirichlet
    U[0]=Ta
    U[N+1]=Tb

    print('La solución es:',U)

    #Se establece el dominio para solución numérica
    XNumerica=np.linspace(0,L,N+2)
    #Se establece el dominio para solución analítica
    XAnalitica = np.linspace(0,L,100)
    #Se calcula la solución analítica
    SolucionA = func.analyticSol(XAnalitica,rho,vel,k,L,Ta,Tb)

    #Salvamos sol
    par['solucion']=U
    par['x']=XNumerica
    hdf5.saveParametros(out_file_name,par)

    #Figuras
    plt.plot(XAnalitica,SolucionA,'-b',label='Solución analítica')
    plt.plot(XNumerica,U,'o--',label='Solución numérica')
    plt.legend()
    plt.grid()

    plt.show()
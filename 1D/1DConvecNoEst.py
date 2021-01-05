# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 11:47:40 2020

@author: Kelly
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import special
from Funciones1D import hdf5 
import sys

def laplaciano1DConveccNOest(N,h,ht,k, rho, vel):
    A=np.zeros((N,N))
    a = (ht * k) / (h**2)
    b = (ht * vel) / (2*h)
    A[0,0]= 2 * a + 1
    A[0,1]= b - a
    A[1,0]= -(a + b)
    for i in range(1,N-1):
        A[i,i] = 2 * a + 1
        A[i+1,i] = -(a + b)
        A[i,i+1] = b - a
    A[N-1,N-2] = -(a + b)
    A[N-1,N-1] = (2 * a + 1)
    return A

def AnalyticSolution(x,t,vel,k):
    den= 2*np.sqrt(k*t)
    solution=0.5*(special.erfc((x-vel*t)/(den)) +
                    np.exp(vel*x)*np.exp(-k)*special.erfc((x+vel*t)/den))
    return solution

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

    par =hdf5.leerParametros(in_file_name,'a','b', 'N', 'Ta', 'Tb', 'k', 'q','ht','Nt','vel','rho','Tol' ) 
    print('Los parámetros empleados son:')

    for key,val in par.items():
            exec(key + '=val')

    L=b-a
    # L = 2.5 # m
    # rho = 1.0 # kg/m^3
    # vel = 1.0 # m/s
    # k = 0.001 # kg / m.s
    # Ta = 1.0 #
    # Tb = 0.0 #
    # N = 300 # Número de incógnitas
    # ht = 0.002
    # Nt = 500
    h = L / (N+1)



    #Solución analítica
    XAnalitica=np.linspace(0,L,100)
    SolA=AnalyticSolution(XAnalitica, ht*Nt, vel, k)

    #RHS

    f = np.zeros(N)


    u = np.zeros(N+2)
    u[0] = Ta
    u[N+1] = Tb

    xnum = np.linspace(0,L,N+2)

    for i in range(1, Nt+1):

        print('Time step = {}'.format(i * ht), sep = '\t')

        A = laplaciano1DConveccNOest(N, h, ht, k, rho, vel)
        f = np.copy(u[1:N+1])
    # Aplicacion de las cond. Dirichlet para el calculo correcto de la sol.
        f[0]   +=  Ta * (rho * vel / (2*h) + k / h**2) * ht
        f[N-1] += -Tb * (rho * vel / (2*h) - k / h**2) * ht

    # Se utiliza un algoritmo del paquete linalg para obtener la solucion del sistema de N x N
        u[1:N+1] = np.linalg.solve(A,f)

        if (i % 100 == 0):
            etiqueta = 'Step = {}'.format(i*ht)
            plt.plot(xnum, u, label=etiqueta)





    plt.plot(XAnalitica,SolA)
    plt.show()
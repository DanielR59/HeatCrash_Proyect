# -*- coding: utf-8 -*-
"""
Created on Thu Dec 31 02:32:02 2020

@author: Kelly
"""

from Funciones1D import funcCondE1D as f
from Funciones1D import hdf5 
import numpy as np 
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

    par =hdf5.leerParametros(in_file_name,'a','b', 'N', 'Ta', 'Tb', 'k', 'q' ) 
    print('Los parámetros empleados son:',par)

    a = par['a']
    b = par['b']
    N = par['N']
    Ta = par['Ta']
    Tb = par['Tb']
    #k = par['k']
    q = par['q']

    x= f.EspaNume(a,b,N)

    #Se calculan k variables 
    k= f.TipoK(x)
    #k =f.promediok(kappa)

    L = b - a
    h = L / (N + 1)
    r = -k / (h**2)
    datos = np.array(10)
    datos = [a,b,N,L,h,Ta,Tb,k,q,r]

        
    k0=k[0]
    kN=k[N+1]


    A1 = f.MatAc3(datos,k)
    #Problema 1
    #Matriz Q
    MatrizQ=f.MatQc3(N,h,q)       # Lado derecho
    #Matriz f 
    MatCondiciones = f.MatDirichletc3(N,h,Ta,Tb,q,k)
    #Matriz b
    Matb = f.Matb2(MatrizQ,MatCondiciones)

    #Parte que soluciona la matriz 
    solucion = f.sol(A1,Matb,datos)
    U = f.u2(solucion,datos)

    #Salvamos parametros
    par['solucion']=U
    par['x']=x
    hdf5.saveParametros(out_file_name,par)
    
    
    
    #GRAFICA
    f.GrafSolc3(x,k,U)







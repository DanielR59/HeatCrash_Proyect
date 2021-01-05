# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 13:15:40 2020

@author: Kelly
"""

import numpy as np
import matplotlib.pyplot as plt
import funcionesCE1D as func
import hdf5CE1D 


#Se leen los parámetros de interés del archivo generado por el programa hdf51D.py desde terminal
par =hdf5CE1D.leerParametros('ENTRADA1','L','rho', 'N', 'kappa', 'vel', 'T0', 'TL' ) 
print(par)

L = par['L']
N = par['N']
h= L/(N+1)
print('El valor de h es:', h)
#kappa = par['kappa']
rho = par['rho']
vel = par['vel']
T0 = par['T0']
TL = par['TL']

#Se establece el dominio para solución numérica

XNumerica=np.linspace(0,L,N+2)


#Se calculan k variables 
k= func.TipoK(XNumerica)
kappa =func.promediok(k)

print('Los valores de K son:', kappa)

#Se establece el dominio para solución analítica

XAnalitica = np.linspace(0,L,100)
SolucionA = func.analyticSolution(XAnalitica,kappa,N,rho,vel,L,T0,TL)

#Se manda a llamar la función para contruir a matriz A
A=func.Laplaciano1DConvEst2(N, h, kappa, rho, vel)
print('La matriz A es:',A)

# Aplicacion de las cond. Dirichlet para el calculo correcto de la sol.
f= np.zeros(N)
f[0] = T0*((rho*vel)/(2*h) + (kappa[0])/(h**2)) 
f[N-1]= -TL*((rho*vel)/(2*h) - (kappa[N-1])/(h**2))
print(f)


# La solucion en el arreglo u, incluye las fronteras
U = np.zeros(N+2)
#obtención de solucion del sistema de N x N

U[1:N+1] = np.linalg.solve(A,f)


#Valores de frontera conocidos debido a las condiciones de Dirichlet

U[0]=T0
U[N+1]=TL

print('La solución es:',U)


#Figuras
func.Graf_Sol_A_N(XNumerica,U,XAnalitica,SolucionA)
func.Graf_k(XNumerica,kappa)

# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 13:15:40 2020

@author: Kelly
"""

import  numpy as np
import matplotlib.pyplot as plt
import funcionesCE1D as func
import hdf5CE1D 

#Se leen los parámetros de interés del archivo generado por el programa hdf51D.py desde terminal
par =hdf5CE1D.leerParametros('ENTRADA1','L','rho', 'N', 'kappa', 'vel', 'T0', 'TL' ) 
print('Los parámetros empleados son:',par)

L = par['L']
N = par['N']
h= L/(N+1)
print('El valor de h es:', h)
kappa = par['kappa']
rho = par['rho']
vel = par['vel']
T0 = par['T0']
TL = par['TL']

#Se manda a llamar la función para contruir a matriz A considerando condiciones
# de frontera tipo Neumann
A=func.Laplaciano1DConvEstNeumann(N, h, kappa, rho, vel)
#Se imprime la matriz A
print('La matriz A es:',A)


# Aplicacion de las cond. Neumann para el calculo correcto de la sol.
f= np.zeros(N+1)
f[0] = T0*((rho*vel)/(2*h) + (kappa)/(h**2)) 
f[N]= -2*h*TL*((rho*vel)/(2*h) - (kappa)/(h**2))
print(f)


# La solucion en el arreglo u, incluye las fronteras
U = np.zeros(N+2)
#obtención de solucion del sistema de N+1 x N+1
U[1:N+2] = np.linalg.solve(A,f)


#Valor de frontera conocido debido a una condición de dirichlet
U[0]=T0
#U[N+1]=TL

print('La solución es:',U)

#Se establece el dominio para solución numérica
XNumerica=np.linspace(0,L,N+2)

print('No se conoce solución analítica con condición de frontera tipo Neumann')


func.Graf_Newmman(XNumerica,U)
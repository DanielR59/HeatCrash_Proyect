

import funcCondE1D as f
import hdf5CondE1D as h
import numpy as np 

par =h.leerParametros('ENTRADA0','a','b', 'N', 'A', 'B', 'k', 'q' ) 
print('Los parámetros empleados son:',par)

a = par['a']
b = par['b']
N = par['N']
A = par['A']
B = par['B']
k = par['k']
q = par['q']

L = b - a
h = L / (N + 1)
r = -k / (h**2)
datos = np.array(10)
datos = [a,b,N,L,h,A,B,k,q,r]


A = f.MatrizA(datos,-2)
x= f.EspaNume(a,b,N)

#Problema 1
#Matriz Q
MatrizQ=f.MatQ(datos)       # Lado derecho
#Matriz f 
MatCondiciones = f.MatDirichlet(datos)
#Matriz b
Matb = f.Matb(MatrizQ,MatCondiciones)

#Parte que soluciona la matriz 
solucion = f.sol(A,Matb,datos)
U = f.u2(solucion,datos)

#Solucion analítica:
xs=f.EspaAna(datos)
sa=f.SolAna1(xs,datos)

#

#GRAFICA
f.GrafSol(x,U,xs,sa)


# -*- coding: utf-8 -*-
"""
Created on Thu Dec 31 02:32:02 2020

@author: Kelly
"""



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
datos = [a,b,N,L,h,A,B,k,q,r]

    
k0=k[0]
kN=k[N+1]


A1 = f.MatAc3(datos,k)
#Problema 1
#Matriz Q
MatrizQ=f.MatQc3(N,h,q)       # Lado derecho
#Matriz f 
MatCondiciones = f.MatDirichletc3(N,h,A,B,q,k)
#Matriz b
Matb = f.Matb2(MatrizQ,MatCondiciones)

#Parte que soluciona la matriz 
solucion = f.sol(A1,Matb,datos)
U = f.u2(solucion,datos)

#GRAFICA
f.GrafSolc3(x,k,U)



# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 13:15:40 2020

@author: Kelly
"""

import  numpy as np
import matplotlib.pyplot as plt
import h5py


def Laplaciano1DConvEst(N, h, kappa, rho, vel):
    """
    

    Parameters
    ----------
    N : Número de incoógnitas
    h : distancia entre nodos
    kappa : conductividad térmica
    rho : Densidad
    vel : componente de velocidad en una dirección
    Returns
    -------
    A : matriz A.

    """
    a = kappa/(h**2)
    b = (rho*vel)/(h)
    A= np.zeros((N,N))
    if vel >0: 
        A[1,0]= -(b + a)
        A[0,0]= b + 2*a 
        A[0,1]= -a
        for i in range(1,N-1):
            A[i+1,i]= -(b + a)
            A[i,i] = b + 2*a
            A[i,i+1]= -a
        A[N-1,N-1] =   b + 2*a  
        A[N-1,N-2] =- ( b + a)
    else:
        A[0,0]= 2*a-b
        A[0,1]= b-a
        A[1,0]= -a
        for i in range(1,N-1):
            A[i,i] = 2*a-b
            A[i+1,i]= -a
            A[i,i+1]= b-a
        A[N-1,N-1] = 2*a-b 
        A[N-1,N-2] =  -a
    return A
    

def analyticSolution(x):
    """
    Función que calcula la solución analítica

    Parameters
    ----------
    x : dominio
    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    return ((np.exp(rho*vel*x / kappa)-1) / (np.exp(rho*vel*L / kappa)-1))*(TL-T0)+T0
    

def f(N):
    f= np.zeros(N)
    if vel>0:
        f[0] = T0*((rho*vel)/(h) + (kappa)/(h**2)) 
        f[N-1]= TL*(kappa/h**2)
    else:
        f[0] = T0*(kappa/h**2) 
        f[N-1]= -TL*((rho*vel)/h - kappa/h**2)
    return f

#PARAMETROS CUALQuIERA
dataset = {'L' :1.0, 'N':20, 'kappa' : 0.1,
        'rho': 1.0, 'vel' : 2.1, 'T0': 1.0, 'TL': 0}

with h5py.File('Datos.hdf5','w') as f:
    for key,value in dataset.items():
        f[key]=value

def leerParametros(archivoHDF5,*parametros,diccionario={}):
    with h5py.File(archivoHDF5,'r') as f:
        for parametro in parametros:
            valor=f.get(parametro) 
            if valor==None:
                print('Parametro :', parametro,' no encontrado')
                continue
                
            if valor.dtype == int:
                diccionario[parametro]=int(np.array(valor))
                
            elif valor.dtype == float:
                diccionario[parametro]=float(np.array(valor))
                
            elif valor.dtype == object:
                diccionario[parametro]=np.array2string(valor) 
        return diccionario
      
#Dependiendo del caso le escribo entre comillas los parametros que quiero 

par =leerParametros('Datos.hdf5','L','rho', 'N', 'kappa', 'vel', 'T0', 'TL' ) 
print(par)

L = par['L']
N = par['N']
h= L/(N+1)

print('El valor de h es:', h)
kappa = par['kappa']
rho = par['rho']
vel = par['vel']
T0 = par['T0']
TL = par['TL']


        
# L= 1.0
# N=20
# h=L /(N+1)
# kappa=0.1
# rho=1.0
# vel=2.1
# T0= 1.0
# TL= 0


       
A=Laplaciano1DConvEst(N, h, kappa, rho, vel)

print(A) 



f=f(N)


#Solucion numerica
U = np.zeros(N+2)

U[1:N+1] = np.linalg.solve(A,f)

#Frontera

U[0]=T0
U[N+1]=TL

print(U)

XNumerica=np.linspace(0,L,N+2)
XAnalitica = np.linspace(0,L,100)
SolucionA = analyticSolution(XAnalitica)
plt.plot(XAnalitica,SolucionA,'-b',label='Solución analítica')
plt.plot(XNumerica,U,'o--',label='Solución numérica')
plt.legend()
plt.grid()

plt.show()
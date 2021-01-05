# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 19:08:46 2020

@author: Kelly
"""
import numpy as np
import matplotlib.pyplot as plt
#### ---------Funciones para convección estacionaria en 1D (Dirichlet)------------------
def Laplaciano1DConvEst(N, h, kappa, rho, vel):
    """   
    Función que realiza el laplaciano para convección estacionara con condiciones
    de frontera tipo Dirichlet
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
    b = (rho*vel)/(2*h)
    A= np.zeros((N,N))
    A[0,0]= 2 * a 
    A[0,1]= b - a
    A[1,0]= - b - a
    for i in range(1,N-1):
        A[i,i] = 2 * a
        A[i+1,i]= - b - a
        A[i,i+1]= b - a
        
    A[N-1,N-1] =  2 * a
    A[N-1,N-2] = - b - a
    return A

def analyticSol(x,rho,vel,kappa,L,T0,TL):
    """
    Función que calcula la solución analítica

    Parameters
    ----------
    x : dominio de la solución analítica
    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    return ((np.exp(rho*vel*x /kappa)-1) / (np.exp(rho*vel*L / kappa)-1))*(TL-T0)+T0
    
#### ---------Funciones para convección estacionaria en 1D (Neumann)------------------

def Laplaciano1DConvEstNeumann(N, h, kappa, rho, vel):
    """
    Función que realiza el laplaciano para convección estacionara con 
    condición de Neumann

    Parameters
    ----------
    N : Número de incógnitas
    h : distancia entre nodos
    kappa : conductividad térmica
    rho : Densidad
    vel : componente de velocidad en una dirección
    Returns
    -------
    A : matriz A.

    """
    a = kappa/(h**2)
    b = (rho*vel)/(2*h)
    A= np.zeros((N+1,N+1))
    A[0,0]= 2 * a 
    A[0,1]= b - a
    A[1,0]= - b - a
    for i in range(1,N):
        A[i,i] = 2 * a
        A[i+1,i]= - b - a
        A[i,i+1]= b - a
    A[N,N-1] = -2*a
    A[N,N] = 2* a
    return A

#### ---------Funciones para convección estacionaria en 1D con k variable (Dirchlet)------------------

def Laplaciano1DConvEst2(N, h, kappa, rho, vel):
    """
    Función que realiza el laplaciano para convección estacionara con condiciones
    de frontera tipo Dirichlet con k variable
    Parameters
    ----------
    N : Número de incógnitas
    h : distancia entre nodos
    kappa : conductividad térmica
    rho : Densidad
    vel : componente de velocidad en una dirección
    Returns
    -------
    A : matriz A.

    """
    #a = kappa/(h**2)
    b = (rho*vel)/(2*h)
    A= np.zeros((N,N))
    A[0,0]= 2 * kappa[0]/(h**2) 
    A[0,1]= b - kappa[0]/(h**2)
    A[1,0]= - b - kappa[1]/(h**2)
    for i in range(1,N-1):
        A[i,i] = 2 * kappa[i]/(h**2)
        A[i+1,i]= - b - kappa[i+1]/(h**2)
        A[i,i+1]= b - kappa[i]/(h**2)
        
    A[N-1,N-1] =  2 * kappa[N-1]/(h**2)
    A[N-1,N-2] = - b - kappa[N-1]/(h**2)
    return A

def analyticSolution(x,kappa,N,rho,vel,L,T0,TL):
    """
    
    Parameters
    ----------
    x : dominio de la solución analítica
    kappa : conducción térmica
    N : Numero de incógnitas
    rho : densidad
    vel : velocidad
    L : Longitud
    T0 : condición de frontera
    TL : Condición de frontera

    Returns
    -------
    sa : resultado de la solución analítica

    """
    for i in range(1,N-1):
        sa=((np.exp(rho*vel*x / kappa[i])-1) / (np.exp(rho*vel*L / kappa[i])-1))*(TL-T0)+T0
    return sa
        
def TipoK(x):
    """
    Se establece el tipo de variabilidad del valor de k

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.

    Returns
    -------
    k : TYPE
        DESCRIPTION.

    """
    opc=int(input('¿El valor de k será...?\n 1. sin(4*pi*x)\n 2. random(x)\n'))
    if opc==1:
        k=np.fabs(np.sin(4*np.pi*x))
    elif opc==2:
        k=np.random.random(size=len(x))
    else:
        print('Opcion no valida')
        exit()
    return k

def promediok(k):
    """
    Obtiene un promedio del valor de k

    Parameters
    ----------
    k : conductividad térmica

    Returns
    -------
    pk : promedio de k

    """
    pk=np.zeros(len(k))
    for i in range(len(pk)-1):
        pk[i]=(k[i]+k[i+1])/2
    return pk

#### ---------Funciones para convección estacionaria en 1D k variable (Neumann)------------------

def Laplaciano1DConvEstNeumann2(N, h, kappa, rho, vel):
    """
    Función que realiza el laplaciano para convección estacionara con 
    condición de Neumann y k variable

    Parameters
    ----------
    N : Número de incógnitas
    h : distancia entre nodos
    kappa : conductividad térmica
    rho : Densidad
    vel : componente de velocidad en una dirección
    Returns
    -------
    A : matriz A.

    """
    
    b = (rho*vel)/(2*h)
    A= np.zeros((N+1,N+1))
    A[0,0]= 2 *  kappa[0]/(h**2)
    A[0,1]= b -  kappa[0]/(h**2)
    A[1,0]= - b -  kappa[1]/(h**2)
    for i in range(1,N):
        A[i,i] = 2 *  kappa[i]/(h**2)
        A[i+1,i]= - b -  kappa[i+1]/(h**2)
        A[i,i+1]= b -  kappa[i]/(h**2)
    A[N,N-1] = -2*  kappa[N]/(h**2)
    A[N,N] = 2*  kappa[N]/(h**2)
    return A

#### ---------Funciones para graficar y guardar imagenes----------------------------
def Graf_Sol_A_N(x,U,xs,sa):
    """
    Genera la grafica que compara solucion analítica a numérica

    Parameters
    ----------
    x : Espacion Numérico.
    U : Solución Numérica.
    xs : Espacio Analítico.
    sa : Solución Analítica. 

    Returns
    -------
    None.

    """
    fig=plt.figure()
    plt.plot(x,U,'-or', label='Solucion numérica')
    plt.plot(xs,sa,'--b',label='Solucion analítica')
    plt.xlabel('Dominio')
    plt.ylabel('Temperatura')
    plt.title('Comparación de soluciones',fontsize=14, color='blue')
    plt.grid()
    plt.legend()
    plt.show()
    Guardar_grafico(fig)

def Graf_k(x,k):
    """
    Genera la grafica del comportamiento de la conductividad

    Parameters
    ----------
    x : Espacio Numérico.
    k : conductividad

    Returns
    -------
    None.

    """
    fig1=plt.figure()
    plt.plot(x,k,'-b')
    plt.title('Comportamiento de la conductividad',fontsize=14, color='blue')
    plt.xlabel('Dominio')
    plt.ylabel('Valor de K')
    plt.grid()
    plt.show()
    Guardar_grafico(fig1)

def Graf_Newmman(x,U):
    """
    Genera la grafica de la solución calculada

    Parameters
    ----------
    x : Espacio Numérico.
    U : Solución numérica

    Returns
    -------
    None.

    """
    fig1=plt.figure()
    plt.plot(x,U,'--or')
    plt.title('Solucion Numérica',fontsize=14, color='blue')
    plt.xlabel('Dominio')
    plt.ylabel('Temperatura')
    plt.grid()
    plt.show()
    Guardar_grafico(fig1)

def Guardar_grafico(Nombre_figura):
    """
    Función que guarda los gráficos realizados en el proceso
    en diferentes formatos

    Parameters
    ----------
    Nombre_figura : figure
        Nombre de la figura a guardar

    Returns
    -------
    None.

    """
    guardar=int(input('\n\tSalvar gráfico?\n1) SI\n2) NO\n: '))
    if guardar==1:
        nombre=input('\nNombre de salida: ')
        print('\n\tSeleccione el formato de salida')
        formato=int(input('\n1) PDF\n2) PNG\n3) JPG\n: '))
        if formato==1:
            Nombre_figura.savefig('{}.pdf'.format(nombre))
            print('\n\tEl gráfico fue salvado con éxito')
        elif formato==2:
            Nombre_figura.savefig('{}.png'.format(nombre))
            print('\n\tEl gráfico fue salvado con éxito')
        elif formato==3:
            Nombre_figura.savefig('{}.jpg'.format(nombre))
            print('\n\tEl gráfico fue salvado con éxito')
        else:
            print('\n\tFormato invalido')
    else:
        print('\n\tEl gráfico no ha sido salvado')

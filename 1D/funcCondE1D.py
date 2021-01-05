# -*- coding: utf-8 -*-
"""
Created on Thu Dec 31 02:02:43 2020

@author: Kelly
"""

import numpy as np
import matplotlib.pyplot as plt


#------- Conduccion estacionaria 1D (Dirichlet)
def MatrizA(d, diag):
    """
    Construye la Matriz A de tamanio NxN

    Parameters
    ----------
    diag : el valor de la diagonal principal que depende de la aproximación
    d : Vector de datos del que se usaran las siguientes posiciones.
    d[2]: Numero de nodos. 

    Returns
    -------
    A : Matriz A del problema
    """
    N=d[2]    
    A = []
    f0 = [diag, 1]
    
    ceros0 = list(map(int, list('0'* (N-2))))
    l0 = f0 + ceros0
    A.append(l0)
    
    cerosm = range(0, N-2)
    
    f1 = [1, diag, 1]
    
    izquierda = []
    for i in cerosm:
      izquierda.append(list(map(int, list('0'*i))))
    derecha = []
    for i in cerosm[::-1]:
      derecha.append(list(map(int, list('0'*i))))
    
    for i in range(len(izquierda)):
      A.append(izquierda[i] + f1 + derecha[i])
    
    linver = l0[::-1]
    A.append(linver)
    print('\n La matriz A es:')
    print(np.matrix(A))
    return A

def EspaNume(a,b,N):
    """
    Hace el espacio Numérico

    Parameters
    ----------
    d : Vector de datos del que se usaran las siguientes posiciones.
    d[0]: Punto inicial a.
    d[1]: Punto final b.
    d[2]: Numero de nodos.

    Returns
    -------
    Espacio Numérico en X

    """
    # a = d[0]
    # b = d[1]
    # N = d[2]    
    x=np.linspace(a,b,N+2)
    return x

def MatQ(d):
    """
    Parameters
    ----------
    d : Vector de datos del que se usaran las siguientes posiciones.
    d[2]: Numero de nodos. 
    d[8]: Valor de q(Escalar). 
    d[9]: Valor de -k / (h**2).
    Returns
    -------
    Devuelve la matriz q

    """
    N=d[2]
    q = np.zeros(N)
    q[1:N-1]=d[8]/d[9]
    print('\n La matriz para Q es:') 
    print(np.array(q))
    return q

def MatDirichlet(d):
    """
    Hace la matriz para condiciones de Dirichlet
    Parameters
    ----------
    d : Vector de datos del que se usaran las siguientes posiciones.
    d[2]: Numero de nodos N. 
    d[5]: Condición Dirichlet en A.
    d[6]: Condición Dirichlet en B.

    Returns
    -------
    Devuelve la matriz f

    """
    N=d[2]
    f = np.zeros(N)
    f[0] = d[5]
    f[N-1] = d[6]
    print('\n La matriz para los valores de condiciones de frontera es:')
    print(np.array(f))
    
    return f

def Matb(Q,f):
    """
    Hace la matriz B
    Parameters
    ----------
    q : Matriz Q
    f : Matriz de condiciones de frontera

    Returns
    -------
    b que es una matriz 

    """
    bmat = Q - f
    print('\n  matriz b es:')
    print(np.array(bmat))
    return bmat

def sol(A,b,d):
    """
    Soluciona el sistema de ecuaciones
    Parameters
    ----------
    A : Matriz A
    b : Matriz b
    d : Vector de datos del que se usaran las siguientes posiciones.
    d[2]: Numero de nodos. 

    Returns
    -------
    u : Matriz solución, sin considerar las temperaturas en las fronteras

    """
    N=d[2]
    u = np.zeros(N+2)
    u[1:N+1] = np.linalg.solve(A,b)
    print('\n La matriz solución sin condiciones de frontera es:')
    print(np.array(u))
    return u

def u2(u,d):
    """
    Parameters
    ----------
    u : Matriz solución sin considerar temperatura en las fronteras
    Ta : Temperatura en la frontera a
    Tb : Temperatura en la frontera b
    N : Número de Nodos
    d : Vector de datos del que se usaran las siguientes posiciones.
    d[2]: Numero de nodos.
    d[5]: Condición Dirichlet en A.
    d[6]: Condición Dirichlet en B.

    Returns
    -------
    u2 : Matriz solución con temperatura en las fronteras

    """
    N=d[2]  
    u2 = []
    u2 = u
    u2[0] = d[5]
    u2[N+1] = d[6]
    print('\n La matriz solución con condiciones de frontera es:')

    print(np.array(u2))
    return u2

def EspaAna(d):
    """
    Hace el espacio Analítico con 100 muestras

    Parameters
    ----------
    d : Vector de datos del que se usaran las siguientes posiciones.
    d[0]: Punto inicial a.
    d[1]: Punto final b.

    Returns
    -------
    Espacio Numérico en X

    """
    a = d[0]
    b = d[1]    
    x=np.linspace(a,b,100)
    return x

def SolAna1(x,d):    
    """
    Funcion que regresa la solución analítica del problema 1. 
    Parameters
    ----------
    x: Espacio generado para evaluar la función
    
    d : Vector de datos del que se usaran las siguientes posiciones.
    En el problema 1:
    d[3]: Longitud del segmento.
    d[5]: Condición Dirichlet en A.
    d[6]: Condición Dirichlet en B.
    d[7]: Valor de k.
    d[8]: Valor de q(Escalar). 

    Returns
    -------
    La solución analítica

    """
    return ((d[6] - d[5])/d[3] + d[8] /(2*d[7]) * (d[3] - x) ) * x + d[5]

def GrafSol(x,U,xs,sa):
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
    plt.title('Comparación de soluciones')
    plt.grid()
    plt.legend(loc='upper left')
    plt.show()
    Guardar_grafico(fig)
    
#------------ K VARIABLE---------------------------

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

def MatAc3(d,pk):
    """
    Parameters
    ----------
    d : Vecto r de Datos
    pk : Vector con promedio de la conductividad no constante

    Returns
    -------
    None.

    """
    N=d[2]    
    A = np.zeros((N,N))
    A[0,0] = -(pk[0]+pk[1]); A[0,1] = pk[1]
    for i in range(1,N-1):
        A[i,i] = -(pk[i]+pk[i+1])
        A[i,i+1] = pk[i+1]
        A[i,i-1] = pk[i]
    A[N-1,N-2] = pk[N-1]; A[N-1,N-1] = -(pk[N]+pk[N-1])
    print('\n La matriz A es:')
    print(np.matrix(A))    
    return A   

def MatQ2(N,k,q):
    """
    Parameters
    ----------
    d : Vector de datos del que se usaran las siguientes posiciones.
    d[2]: Numero de nodos. 
    d[8]: Valor de q(Escalar). 
    d[9]: Valor de -k / (h**2).
    Returns
    -------
    Devuelve la matriz q

    """
    for i in range(1,N-1):
        q1 = np.zeros(N)
        q1[1:N-1]=q/k[i]
    print('\n La matriz para Q es:') 
    print(np.array(q))
    return q1

def GrafSol2(x,U,k):
    """
    Genera la grafica de la solucion analítica 
    y comportamiento de la conductividad

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
    fig1=plt.figure()
    plt.plot(x,U,'-or', label='Solucion numérica')
    plt.title('Solucion')
    plt.xlabel('Dominio')
    plt.ylabel('Temperatura')
    plt.legend(loc='upper left')
    plt.grid()
    plt.show()
    Guardar_grafico(fig1)
    
    fig2=plt.figure(2)
    plt.plot(x,k,'-b',label='K variable')
    plt.title('Solucion')
    plt.xlabel('Dominio')
    plt.ylabel('Valor de K')
    plt.legend(loc='upper left')
    plt.grid()
    plt.show()
    Guardar_grafico(fig2)
    
def MatQc3(N,h,q):
    """
    Parameters
    ----------
    d : Vector de datos del que se usaran las siguientes posiciones.
    d[2]: Numero de nodos.
    d[4]: Tamanio de la malla h.
    d[7]: Valor de q(Escalar). 
    
    Returns
    -------
    Devuelve la matriz q

    """
   # N=d[2]
    q1 = np.zeros(N)
    q1[1:N-1]=q*h*h
    print('\n La matriz para Q es:') 
    print(np.array(q))
    return q1

def MatDirichletc3(N,h,A,B,q,k):
    """
    Hace la matriz para condiciones de Dirichlet
    Parameters
    ----------
    d : Vector de datos del que se usaran las siguientes posiciones.
    d[2]: Numero de nodos N. 
    d[4]: Tamanio de la malla h.
    d[5]: Condición Dirichlet en A.
    d[6]: Condición Dirichlet en B.
    d[8]: Valor de q(Escalar).
    
    Returns
    -------
    Devuelve la matriz f

    """
    #N=d[2]
    f = np.zeros(N)
    f[0] = (h**2)*q-k[0]*A
    f[N-1] = h*h*q-k[N]*B
    print('\n La matriz para los valores de condiciones de frontera es:')
    print(np.array(f))
    
    return f

def Matb2(Q,f):
    """
    Hace la matriz B
    Parameters
    ----------
    q : Matriz Q
    f : Matriz de condiciones de frontera

    Returns
    -------
    b que es una matriz 

    """
    bmat = Q + f
    print('\n  matriz b es:')
    print(np.array(bmat))
    return bmat

def GrafSolc3(x,k,u):
    """
    
    Parameters
    ----------
    x : Dominio numérico
    k : Valor de la conductividad no constante
    u : Solución numérica

    Returns
    -------
    None.

    """
    fig, ax=plt.subplots(2,figsize=(12,8))
    fig.suptitle('Caso con conductividad no constante', fontsize =18)
    
    ax[0].plot(x,k,'r')
    ax[0].set_title('Comportamiento de K',color='blue')
    ax[0].set(ylabel='K')
    ax[0].grid(True)
    
    ax[1].plot(x,u,'-bo')
    ax[1].set_title('Solucion numerica',color='blue')
    ax[1].set(xlabel='Dominio', ylabel='Solución')
    ax[1].grid(True)
    
    plt.show()
    
    Guardar_grafico(fig)
    
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

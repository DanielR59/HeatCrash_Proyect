# PROGRAMA QUE CONTIENE LAS FUNCIONES PARA EL CÁLCULO DE LA SOLUCIÓN
# DE LA ECUACIÓN DE CALOR EN ESTADO NO ESTACIONARIO Y CONSIDERANDO EL 
# TÉRMINO CONVECTIVO

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import special

def lectura_variables(sel,Titulo):
    """
    Esta función permite seleccionar al usuario la manera 
    en la que se ingresarán los datos

    Parameters
    ----------
    sel : Integer
        Valor de selección ingresado por el usuario

    Returns
    -------
    a : float
        Valor al inicio del dominio.
    b : float
        Valor al final del dominio.
    N : Integer
        Número de nodos
    Ta : float
        Temperatura en la frontera del inicio.
    Tb : float
        Temperatura en la frontera del final.
    S : float
        Valor de la fuente
    k : float
        Valor de la conductividad térmica
    ht : float
        Tamaño de paso en el tiempo
    Tmax: float
        Tiempo máximo hasta el cual calculamos la solución
    v : floar
        Velocidad del fluido en movimiento
    Tol : float
        Tolerancia del algoritmo
        

    """
    if sel == 1:
        val= np.loadtxt(Titulo)
        a = val[0]
        b = val[1]
        N = int(val[2])
        Ta = val[3]
        Tb = val[4]
        k = val[5]
        S = val[6]
        ht = val[7]
        Tmax = val[8]
        v = val[9]
        Tol = val[10]
       
        return a,b,N,Ta,Tb,k,S,ht,Tmax,v,Tol
    else:    
        # Datos de entrada
        a = float(input("\tIngrese el comienzo de la barra.\na= "))
        b = float(input("\tIngrese el fin de la barra.\nb= "))
        N = int(input("\tIngresa el número de nodos que desea.\nN="))
        Ta = float(input("\tIngrese la temperaruta al inicio.\nTa="))
        Tb = float(input("\tIngrese la temperaruta al final.\nTb="))
        k = float(input("\tIngrese la conductividad térmica.\nk="))
        S = float(input("\tIngrese las fuentes o sumideros.\nS="))
        ht = float(input("\tIngrese el paso del tiempo.\nht="))
        Tmax= float(input("\tIngrese el tiempo maximo.\nTmax="))
        v = float(input("\tIngrese la velocidad.\nv="))
        Tol = float(input("\tIngrese la Tolerancia.\nTol="))
        return a,b,N,Ta,Tb,k,S, ht, Tmax,v,Tol

def constantes(a,b,N,ht,Tmax,K,v):
    """
    Esta función calcula las constantes que requiere el programa

    Parameters
    ----------
    a : float
        Valor al inicio del dominio.
    b : float
        Valor al final del dominio.
    N : Integer
        Número de nodos.
    ht : float
         Paso en el tiempo
    k : float
        Valor de la conductividad térmica
    Tmax: float
        Tiempo máximo hasta el cual calculamos la solución
    v : float
        Velocidad del fluido en movimiento
         
    Returns
    -------
    h : float
        Distancia entre cada nodo
    x : array
        Vector con el cual se graficará
    lar : float
        Distancia total del dominio
    Nt: Int
        Número de nodos en el tiempo
    r : float
    p : float

    """
#    v = float(v)
    h = (b-a)/(N+1)
    x = np.linspace(a,b,N+2)
    lar = b-a
    Nt = int(Tmax / ht)
    r = ht * K / (h*h)
    p = ht *  v / (2*h)
    return h,x,lar,Nt,r, p

def dominio(a,b,N):
    """
    Esta función genera la malla del dominio.

    Parameters
    ----------
    a : Float
        Incio del dominio.
    b : Float
        Fin del dominio.
    N : Int
        Numero de nodos.

    Returns
    -------
    x : array
        Vector de la malla del dominio.

    """
    x = np.linspace(a,b,N+2)
    return x

def vector_aux(Ta,Tb,N,q):
    """
    Esta función genera un vector el cual se utiliza en la resolución del 
    problema.

    Parameters
    ----------
    Ta : float
        Temperatura en la frontera del inicio.
    Tb : float
        Temperatura en la frontera del final.
    N : Integer
        Número de nodos.
    q : array
        Vector con la información de las fuentes
    Returns
    -------
    b + q: array
        Vector con las condiciones a la frontera

    """
    b = np.zeros(N+2)
    b[0] = Ta
    b[-1] = Tb
    return b + q


def matriz_Diagonal_Conv_Noest(N,r,p):
    """
    Esta función genera la matriz diagonal necesaria en la resolución del 
    problema.    

    Parameters
    ----------
    N : Integer
        Número de nodos.
    r : Float
        Valor en la diagonal
    p : Float
        Valor en la diagonal

    Returns
    -------
    A : np array
        Matriz diagonal

    """
    A = np.zeros((N,N))
    
    A[0,0] = 1 + (2 * r); 
    A[0,1] = p-r
    for i in range(1,N-1):
        A[i,i] = 1 + (2 * r) 
        A[i,i+1] = p-r
        A[i,i-1] = -(p+r)
    A[N-1,N-2] = -(p+r); 
    A[N-1,N-1] = 1 + (2 * r)
    return A


def graficas(xa,ua,titulo):
    """
    Esta función genera las gráficas de la solución analítica y numérica

    Parameters
    ----------
    xa : float
        Vector a graficar en X
    ua : float
        Vector a graficar en Y
    titulo : str
        Cadena con el título de la gráfica

    Returns
    -------
    Grafica de la solución.
    """
    plt.plot(xa,ua, 'k-', lw=2.5, label='Solución Analítica')
    plt.xlabel('$x$')
    plt.ylabel('$u(x)$')
    plt.title(titulo)
    plt.grid()
    plt.legend()
    plt.show()
    
def grafica_error(Error):
    """
    Esta función grafica el error calculado.

    Parameters
    ----------
    Error : float
        Vector que contiene el error de la solución.

    Returns
    -------
    None.

    """
    plt.plot(Error)
    plt.title('Error en la solución')
    plt.yscale('log')
    plt.xlabel('$n$')
    plt.grid()
    plt.ylabel('$RMS$')
    plt.show()
    

    
def sol_Analitica(a,b,K,Tmax,N,v):  
    """
    Esta función genera un vector que contiene la solución exacta al problema.

    Parameters
    ----------
    a : float
        Valor al inicio del dominio.
    b : float
        Valor al final del dominio.
    N : Integer
        Número de nodos.
    k : float
        Valor de la conductividad térmica
    Tmax: float
        Tiempo máximo hasta el cual calculamos la solución
    v : floar
        Velocidad del fluido en movimiento
         
    Returns
    -------
    xa : float.
        Malla del dominio sin incluir las fronteras.
    ua : float.
        Solución exacta al problema.

    """
    xa = np.linspace(a, b, N)
    
    divisor = 2 * np.sqrt(K * Tmax)
    ua = 0.5 * (special.erfc((xa - v * Tmax)/ divisor) + 
                np.exp(v * xa) * np.exp(-K) * special.erfc((xa + v * Tmax)/divisor))
    return (xa, ua)
    
def escritura(u,u_exa):
    """
    Esta función genera un archivo .txt de dos columnas con los datos de
    las soluciones numérica y analítica
    Parameters
    ----------
    u : float
        Vector que contiene la solución numérica del problema
    u_exa : float
        Vector que contiene la solución analítica del problema

    Returns
    -------
    Archivo.
    """
    
    serie1 = pd.Series(u)
    serie2 = pd.Series(u_exa)
    tabla = pd.DataFrame(serie1,columns = ['Solución analítica'])
    tabla['Solución numérica'] = serie2
    np.savetxt('Solución1.txt',tabla,fmt='%f', header = 'Sol. Analítica    Sol. Exacta')
    
    

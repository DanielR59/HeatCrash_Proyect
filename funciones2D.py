# -*- coding: utf-8 -*-
"""



@author: daniel
"""
import numpy as np

# =============================================================================
# Condiciones de frontera
# =============================================================================
def boundary_cond_dirichtlet(matriz,Tx1,Tx2,Ty1,Ty2):
    """
    Funcion que agrega las condiciones de frontera tipo Dirichlet a la matriz
                           Tx2
      0,0                 frontera A
       o-----o-----o-----o-----o-----o-----o-----o
    Ty1|     |     |     |     |     |     |     |Ty2  f
       o-----o-----o-----o-----o-----o-----o-----o     r
       |     |     |     |     |     |     |     |     o
     C o-----o-----o-----o-----o-----o-----o-----o     n
     C |     |     |     |     |     |     |     | Ny  t
     C o-----o-----o-----o-----o-----o-----o-----o     e
       |     |     |     |     |       |   |     |
       o-----o-----o-----o-----o-----o-----o-----o    D
       |     |     |     |     |     |     |     |
       o-----o-----o-----o-----o-----o-----o-----o
                           Nx
                           Tx1
                       frontera B

    Parameters
    ----------
    matriz : numpy array
        DESCRIPTION.
    Tx1 : float
        DESCRIPTION.
    Tx2 : float
        DESCRIPTION.
    Ty1 : float
        DESCRIPTION.
    Ty2 : float
        DESCRIPTION.

    Returns
    -------
    matriz : TYPE
        DESCRIPTION.

    """
    matriz[-1,:] = Tx2
    matriz[:,0] = Ty1
    matriz[:,-1] = Ty2
    matriz[0,:] = Tx1
    return matriz



# =============================================================================
# Funciones para casos Estacionarios
# =============================================================================

def left_matrix_2d(Nx,Ny):
    """
    Funcion que genera la matriz cuadrada izquierda de dimensiones Nx*Ny para
    resolver el sistema del caso estacionario en 2D

    Parameters
    ----------
    Nx : int
        DESCRIPTION.
    Ny : int
        DESCRIPTION.

    Returns
    -------
    matriz_aux1 : numpy array
        DESCRIPTION.

    """
    multiples_auxiliares=np.array(range(Nx,Nx*Ny,Nx))#No preguntes solo gozalo
    multiples_auxiliares-=1#Porque python empieza en 0

    aux=-4 #valor de la diagonal principal

    matriz_aux1=np.zeros([Nx*Ny,Nx*Ny])#creo matriz Nx*Ny con diagonal principal
    for i in range(0,Nx*Ny):
        matriz_aux1[i,i]=aux
    for i in range(0,Nx*Ny-1): #arreglo que llena los valores proximos a la diagonal principal

        if i not in multiples_auxiliares:
            matriz_aux1[i+1,i]=1
            matriz_aux1[i,i+1]=1

    for i in range(0,Nx*Ny-Nx):#Se llenan los valores de la matriz
        matriz_aux1[Nx+i,i]=1

    for i in range(0,Nx*Ny-Nx):
        matriz_aux1[i,Nx+i]=1
    #
    # matriz_aux1/=(hx*hy)
    # print(matriz_aux1)
    return matriz_aux1


# =============================================================================
# Funciones para casos NO Estacionarios
# =============================================================================
def iterationTime2D(u,q,hx,hy,ht,k):
    """
    Funcion que da una iteracion en tiempo para resolver los casos de Conduccion
    de Calor no estacionario mediante el método de Euler implicito

    Parameters
    ----------
    u : numpy array
        DESCRIPTION.
    q : numpy array
        DESCRIPTION.
    hx : float
        DESCRIPTION.
    hy : float
        DESCRIPTION.
    ht : float
        DESCRIPTION.
    k : float
        DESCRIPTION.

    Returns
    -------
    u_updated : numpy array
        DESCRIPTION.
    error : float
        DESCRIPTION.

    """


    u_updated=u.copy()

    u_updated[1:-1,1:-1] = u[1:-1,1:-1]+\
        (k * ht / hy**2)*(u[1:-1, 2:] - 2 * u[1:-1, 1:-1] + u[1:-1, :-2])+\
            (k * ht / hx**2)* (u[2:,1: -1] - 2 * u[1:-1, 1:-1] + u[:-2, 1:-1])+ht*q[1:-1,1:-1]
    error=np.linalg.norm(u_updated-u,2)
    return u_updated,error

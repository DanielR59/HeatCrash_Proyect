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

def iterationCond2D(u,q,hx,hy,kx,ky):
    """
    iterationCond2D [summary]

    [extended_summary]

    Parameters
    ----------
    u : [type]
        [description]
    q : [type]
        [description]
    hx : [type]
        [description]
    hy : [type]
        [description]
    kx : [type]
        [description]
    ky : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    B = kx/hx**2
    C = kx/hx**2
    D = ky/hy**2
    E = ky/hy**2
    A= -(B+C+D+E)
    u_updated = u.copy()
    u_updated[1:-1,1:-1] = (q[1:-1,1:-1]-B*u[:-2,1:-1]-C*u[2:,1:-1]-D*u[1:-1,:-2]-E*u[1:-1,2:])/A
    # u_updated[1:-1,1:-1] =(hy**2*(u[:-2,1:-1] + u[2:,1:-1]) \
        # + hx**2*(u[1:-1,:-2] + u[1:-1,2:])- hy**2*hx**2*q[1:-1,1:-1])/(2*(hx**2+hy**2))
    
    error=np.linalg.norm(u_updated-u,2)
    return u_updated,error
def iterationConv2D(u,q,alpha_x,alpha_y,kappa_x,kappa_y,hx,hy):
    """
    iterationConv2D [summary]

    [extended_summary]

    Parameters
    ----------
    u : [type]
        [description]
    q : [type]
        [description]
    alpha_x : [type]
        [description]
    alpha_y : [type]
        [description]
    kappa_x : [type]
        [description]
    kappa_y : [type]
        [description]
    hx : [type]
        [description]
    hy : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    aux_diagp=2*(kappa_x/hx**2+kappa_y/hy**2) #valor de la diagonal principal

    aux_1i = alpha_x/(2*hx)
    aux_2i = kappa_x/hx**2
    aux_1j = alpha_y/(2*hy)
    aux_2j = kappa_y/hy**2

    u_updated=u.copy()

    u_updated[1:-1,1:-1] = (u[2:,1:-1]*(aux_2i-aux_1i)+u[:-2,1:-1]*(aux_2i+aux_1i)+\
        u[1:-1,2:]*(aux_2j-aux_1j)+u[1:-1,:-2]*(aux_2j+aux_1j)+q[1:-1,1:-1])/(2*aux_2i+2*aux_2j)


    error=np.linalg.norm(u_updated-u,2)
    return u_updated,error

# =============================================================================
# Funciones para casos NO Estacionarios
# =============================================================================
def iterationTime2D(u,q,hx,hy,ht,kappa_x,kappa_y):
    """
    Funcion que da una iteracion en tiempo para resolver los casos de Conduccion
    de Calor no estacionario mediante el método de Euler 

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
        (kappa_y * ht / hy**2)*(u[1:-1, 2:] - 2 * u[1:-1, 1:-1] + u[1:-1, :-2])+\
            (kappa_x * ht / hx**2)* (u[2:,1: -1] - 2 * u[1:-1, 1:-1] + u[:-2, 1:-1])+ht*q[1:-1,1:-1]
    error=np.linalg.norm(u_updated-u,2)
    return u_updated,error


def iterationtimeConvNoEst2D(u,q,alpha_x,alpha_y,kappa_x,kappa_y,hx,hy,ht):
    """
    iterationtimeConvNoEst2D [summary]

    [extended_summary]

    Parameters
    ----------
    u : [type]
        [description]
    q : [type]
        [description]
    alpha_x : [type]
        [description]
    alpha_y : [type]
        [description]
    kappa_x : [type]
        [description]
    kappa_y : [type]
        [description]
    hx : [type]
        [description]
    hy : [type]
        [description]
    ht : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    
    u_updated=u.copy()
    u_updated[1:-1,1:-1] = u[1:-1,1:-1]+\
            (kappa_y * ht / hy**2)*(u[1:-1, 2:] - 2 * u[1:-1, 1:-1] + u[1:-1, :-2])+\
            (kappa_x * ht / hx**2)* (u[2:,1: -1] - 2 * u[1:-1, 1:-1] + u[:-2, 1:-1])+\
            (alpha_x*ht/(2*hx))*(u[1:-1, 2:]-u[1:-1, :-2])+ht*q[1:-1,1:-1]+(alpha_y*ht/(2*hy))*(u[1:-1, 2:]-u[1:-1, :-2])
    

    
    error=np.linalg.norm(u_updated-u,2)
    return u_updated,error
# -*- coding: utf-8 -*-
"""



@author: daniel
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

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

# =============================================================================
# Función para crear gráficos
# =============================================================================
def GrafError(error,color):
    """
    Función para graficar el error

    Parameters
    ----------
    error : numpy array
        Arreglo con los valores de los errores
    color : char
        Color de la linea de trazado.

    Returns
    -------
    None.

    """
    err=plt.figure(figsize=(5,5))
    plt.plot(error, color)
    plt.semilogy()
    plt.xlabel('Pasos de tiempo')
    plt.ylabel('Error')
    plt.title('Gráfica de tendencia del error')
    plt.grid()
    plt.show()
    Guardar_grafico(err)
    
def GrafContornos(x,y,z,niveles,transparencia,mapaColor):
    """
    Función que realiza un gráfico de contornos

    Parameters
    ----------
    x : numpy array
        Datos del dominio en x.
    y : numpy array
        Datos del dominio en y.
    z : numpy array
        Matriz con datos de la malla.
    niveles : int
        Numero de contornos a usar
    transparencia : float
        Transparencia del mapa, valores entre 0 y 1
    mapaColor :string
        Tipo de mapa de color a usar

    Returns
    -------
    None.

    """
    contorno=plt.figure(figsize=(8,5))  
    c=plt.contourf(x, y, z,niveles, alpha=transparencia,cmap=mapaColor)
    cbar=contorno.colorbar(c, shrink=1.0)
    cbar.set_label('Temperatura [ºC]')
    plt.xlabel('Distancia x [m]')
    plt.ylabel('Distancia y [m]')
    plt.title('Distribución de la temperatura (mapa de contornos)')
    plt.show()
    Guardar_grafico(contorno)
       
def Graf3D(x,y,z,mapaColor):
    """
    Función que realiza un gráfico en 3D

    Parameters
    ----------
    x : numpy array
        Datos del dominio en x.
    y : numpy array
        Datos del dominio en y.
    z : numpy array
        Matriz con datos del dominio en z
    mapaColor :string
        Tipo de mapa de color a usar

    Returns
    -------
    None.

    """
    surf=plt.figure(figsize=(5,4)) 
    ax=surf.gca(projection='3d')
    s=ax.plot_surface(x, y, z, cmap=mapaColor)
    cbar=surf.colorbar(s, shrink=0.5)
    cbar.set_label('Temperatura [ºC]')
    plt.xlabel('Distancia x [m]')
    plt.ylabel('Distancia y [m]')
    plt.title('Distribución de la temperatura (superficie 3D)')
    plt.show()
    Guardar_grafico(surf)
    
def Animacion_Contorno(x,y,z,Nt):
    """
    Funcion de animación de mapas de contorno

    Parameters
    ----------
    x : Array
        Dominio x
    y : Array
        Dominio y
    z : Array
        Dominio z
    Nt : float
        Pasos de tiempo

    Returns
    -------
    None.

    """
    contorno=plt.figure(figsize=(8,5))  
    c=plt.contourf(x, y, z[0],8, alpha=0.75,cmap='inferno')
    cbar=contorno.colorbar(c, shrink=1.0)
    cbar.set_label('Temperatura [ºC]')
    plt.xlabel('Distancia x [m]')
    plt.ylabel('Distancia y [m]')
    plt.title('Distribución de la temperatura (mapa de contornos)')
    
    
    
    def cambio(i):
        plt.contourf(x, y, z[i+1], 8, alpha=0.75,cmap='inferno')
            
    anim=animation.FuncAnimation(contorno,            # La figura
                     cambio, # la función que cambia los datos
                     interval=1,     # Intervalo entre cuadros en milisegundos
                     frames=Nt-1,   # Cuadros
                     repeat=False)   # Permite poner la animación en un ciclo
    plt.show()
    answer=input('Quieres guardar la animación?  [y/n]\n')
    if answer=='y':
        Nom_Anim=input('Nombre de la animación : ')
        anim.save('{}.gif'.format(Nom_Anim), writer='imagemagick', fps=10)
    elif answer =='n':
        print('\n\tNo guardado')
    else:
        print('\n\tNo guardado')
    
   
def Animacion_Superficie(x,y,z,Nt):
    """
    Funcion de animación de mapas 3D
    Parameters
    ----------
    x : Array
        Dominio x
    y : Array
        Dominio y
    z : Array
        Dominio z
    Nt : float
        Pasos de tiempo

    Returns
    -------
    None.

    """
    surf=plt.figure(figsize=(5,4)) 
    ax=surf.gca(projection='3d')
    s=ax.plot_surface(x, y, z[0], cmap='inferno')
    cbar=surf.colorbar(s, shrink=0.5)
    cbar.set_label('Temperatura [ºC]')
    plt.xlabel('Distancia x [m]')
    plt.ylabel('Distancia y [m]')
    plt.title('Distribución de la temperatura (superficie 3D)')
    
    def cambio(i):
        ax.plot_surface(x, y, z[i+1], cmap='inferno')
            
    anim=animation.FuncAnimation(surf,            # La figura
                     cambio, # la función que cambia los datos
                     interval=1,     # Intervalo entre cuadros en milisegundos
                     frames=Nt-1,   # Cuadros
                     repeat=False)   # Permite poner la animación en un ciclo
    plt.show()
    answer=input('Quieres guardar la animación?  [y/n]\n')
    if answer=='y':
        Nom_Anim=input('Nombre de la animación : ')
        anim.save('{}.gif'.format(Nom_Anim), writer='imagemagick', fps=10)
    elif answer =='n':
        print('\n\tNo guardado')
    else:
        print('\n\tNo guardado')
    
    
    
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
    

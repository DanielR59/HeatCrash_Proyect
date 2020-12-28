import hdf5
import numpy as np
import sys
import time
import matplotlib.pyplot as plt
from funciones2D import  boundary_cond_dirichtlet,iterationCond2D

if __name__ == "__main__":
    #                       Tx2
    #  0,0                 frontera A
    #   o-----o-----o-----o-----o-----o-----o-----o
    #   |     |     |     |     |     |     |     |     f
    #   o-----o-----o-----o-----o-----o-----o-----o     r
    #   |     |     |     |     |     |     |     |     o
    # C o-----o-----o-----o-----o-----o-----o-----o     n
    # C |     |     |     |     |     |     |     | Ny  t
    # C o-----o-----o-----o-----o-----o-----o-----o     e
    #   |     |     |     |     |       |   |     |
    #   o-----o-----o-----o-----o-----o-----o-----o    D
    #   |     |     |     |     |     |     |     |
    #   o-----o-----o-----o-----o-----o-----o-----o
    #                       Nx
    #                       Tx1

    #                   frontera B
    try:
        in_file_name = sys.argv[1]; out_file_name = sys.argv[2]
    except:
        mensaje = """ Error: La ejecucion de este programa requiere de 2 argumentos.
        Ejecucion correcta: python {} entrada salida
        donde "entrada" es el nombre de un archivo que contiene los
        datos del problema :  este se puede generar con el programa hdf5.py.
        El nombre "salida" se usa para almacenar la solucion del problema.

        Por ejemplo: python {} INPUT_03 SALIDA""".format(__file__,__file__)

        print(mensaje)
        sys.exit(1)

    Datos=hdf5.leerParametros(in_file_name,'ax','ay','bx','by','Nx','Ny','Tx1','Tx2','Ty1','Ty2','Tolerancia')

    for key,val in Datos.items():
        print(key,'=',val)
        print('-'*10)
        exec(key + '=val')


    hx = (bx-ax)/(Nx+1)
    hy = (by-ay)/(Ny+1)
    x = np.linspace(ax,bx,Nx+2)
    y = np.linspace(ay,by,Ny+2)
    xg, yg = np.meshgrid(x,y)

    print('hx = ',hx)
    print('hy = ',hy)

    u = np.zeros((Ny+2, Nx+2))
    u= boundary_cond_dirichtlet(u,Tx1,Tx2,Ty1,Ty2)
    f = np.ones_like(u)*0 # RHS
    
    for i in range(20000):
        u,error=iterationCond2D(u,f,hx,hy)
        if error < Tolerancia:
            break

    
    Datos['xg']=xg
    Datos['yg']=yg
    Datos['solucion']=u

    hdf5.saveParametros(out_file_name,Datos)


    f1=plt.figure()
    c= plt.contourf(xg, yg, u, 8, alpha=.75,cmap='inferno')
    f1.colorbar(c, shrink=1.0)

    f2 = plt.figure()
    ax = f2.gca(projection='3d')
    s = ax.plot_surface(xg, yg, u, cmap='inferno')
    f2.colorbar(s, shrink=0.5)

    # plt.show()
    plt.show()
    
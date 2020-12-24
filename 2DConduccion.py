import hdf5
import numpy as np
import sys
import time
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix



def left_matrix_2D(Nx,Ny):
        
    multiples_auxiliares=np.array([i for i in range(Nx,Nx*Ny,Nx)]) #No preguntes solo gozalo
    multiples_auxiliares-=1 #Porque python empieza en 0 

    aux=-4 #valor de la diagonal principal
    
    matriz_aux1=np.zeros([Nx*Ny,Nx*Ny]) #creo matriz Nx*Ny con diagonal principal
    for i in range(0,Nx*Ny):
        matriz_aux1[i,i]=aux
    for i in range(0,Nx*Ny-1): #arreglo que llena los valores proximos a la diagonal principal
        
        if (i not in multiples_auxiliares):
            matriz_aux1[i+1,i]=1
            matriz_aux1[i,i+1]=1

    for i in range(0,Nx*Ny-Nx): #Se llenan los valores de la matriz 
        matriz_aux1[Nx+i,i]=1
        
    for i in range(0,Nx*Ny-Nx):
        matriz_aux1[i,Nx+i]=1
    # 
    # matriz_aux1/=(hx*hy)
    # print(matriz_aux1)
    return matriz_aux1





if __name__ == "__main__":
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

    #                   frontera B
    try:
        in_file_name = sys.argv[1]; out_file_name = sys.argv[2]
    except:
        mensaje = """ Error: La ejecucion de este programa requiere de 2 argumentos.
        Ejecucion correcta: python 2D_Poisson_01.py entrada salida
        donde "entrada" es el nombre de un archivo que contiene los
        datos del problema :  este se puede generar con el programa hdf5.py.
        El nombre "salida" se usa para almacenar la solucion del problema.

        Por ejemplo: python 2D_Poisson_01.py INPUT_03 SALIDA"""

        print(mensaje)
        sys.exit(1)

    
    
    Datos=hdf5.leerParametros(in_file_name,'ax','ay','bx','by','Nx','Ny','Tx1','Tx2','Ty1','Ty2')
    
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

    f = np.ones((Ny,Nx))*0 # RHS
    # for i in range(1,Nx):
    #     f[:,i]=50*np.sin(x[i])
    f[20,20]=200
    print(f[20,20])
    f[0,:]-=Tx1
    f[-1,:]-=Tx2
    f[:,0]-=Ty1
    f[:,-1]-=Ty2
    
    
    # A=csr_matrix(Laplaciano2D(Nx,Ny,-4)) #Convertimos a una matriz sparse
    A=csr_matrix(left_matrix_2D(Nx,Ny))
    # A=csr_matrix(left_matrix_2D(hx,hy,Nx,Ny))
    # print(A.todense())
    # print(C.todense())
    # print((A.todense()==B.todense()).all())
    # print((B.todense()==C.todense()).all())

    u = np.zeros((Ny+2, Nx+2))
    u[-1,:   ] = Tx2 
    u[:   ,0   ] = Ty1 
    u[:   ,-1] = Ty2
    u[0   ,:   ] = Tx1 
    
    u_sistema = np.zeros([Ny*Nx])
    # print(A.shape,f.flatten().shape)
    print('-'*80)
    t1_start = time.process_time()
    u_sistema=spsolve(A,f.flatten())
    t1_stop = time.process_time()
    print('Tiempo para resolver el sistema = {:0.6f} s'.format(t1_stop-t1_start))
    print('-'*80)
    u_sistema.shape = (Ny, Nx) 

    
    u[1:Ny+1,1:Nx+1] = u_sistema
    
    Datos['xgrid']=xg
    Datos['ygrid']=yg
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
    
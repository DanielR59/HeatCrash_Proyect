import hdf5
import numpy as np
import sys
import time
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from funciones2D import iterationTime2D, boundaryCondDirichtlet

if __name__ == "__main__":
    
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

    
    #leemos los parametros
    Datos=hdf5.leerParametros(in_file_name,'ax','ay','bx','by','Nx','Ny',\
        'Tx1','Tx2','Ty1','Ty2','ht','Tolerancia','Tmax','ht','Tini')
    
    #imprimimos los parametros y los evaluamos
    print('Parametros')
    print('-'*10)
    for key,val in Datos.items():
        print(key,'=',val)
        print('-'*10)
        exec(key + '=val')
        

    hx = (bx-ax)/(Nx+1)
    hy = (by-ay)/(Ny+1)
    
    

    #Checamos la estabilidad
    k=1


    if k*ht/hx**2+k*ht/hy**2<0.5:
        pass
    else:

        print('Metodo inestable se cambia el paso en tiempo para cumplir estabilidad')
        print('ht : ',ht,'--->',min(hx**2,hy**2)/(2*k))
        ht=0.5/k*(1/hx**2+1/hy**2)**-1
        
        


    Nt=int(Tmax/ht)
    print('Nt = ',Nt)
    print('hx = ',hx)
    print('hy = ',hy)
    #Generamos dominio
    x = np.linspace(ax,bx,Nx+2)
    y = np.linspace(ay,by,Ny+2)
    xg, yg = np.meshgrid(x,y)

    #Generamos matriz

    u=np.ones([Ny+2,Nx+2])*Tini #Aplicamos la condicion inicial a la matriz
    #Aplicamos las condiciones de frontera
    u= boundaryCondDirichtlet(u,Tx1,Tx2,Ty1,Ty2)


    q=np.ones_like(u)*0
    # q[5,5]=200500

    errores=[]
    solucion=np.empty([50000,Ny+2,Nx+2],dtype=np.float16)
    for n in range(Nt+1):
        solucion[n,:,:]=u
        u,error=iterationTime2D(u,q,hx,hy,ht,k)
        errores.append(error)
        if error < Tolerancia:
            print('Iteracion terminada con ',n,' pasos')
            break
    
    solucion=solucion[:n+1,:,:]
    
    answer=input('Quieres guardar la solucion para generar una animación?  [y/n]\n')
    if answer =='y':
        Datos['solucion_animada']=solucion
    elif answer =='n':
        pass
    else:
        pass
    Datos['error']=errores
    Datos['xg']=xg
    Datos['yg']=yg
    Datos['solucion']=u
    Datos['Nt']=n
    Datos['ht']=ht
    hdf5.saveParametros(out_file_name,Datos)

    fig = plt.figure(figsize=(11, 7), dpi=100)
    plt.contourf(xg, yg, u)
    plt.colorbar()
    plt.contour(xg, yg, u)
    plt.xlabel('X')
    plt.ylabel('Y')
    
    plt.figure()
    plt.plot(errores)
    plt.semilogy()
    
    
    plt.show()
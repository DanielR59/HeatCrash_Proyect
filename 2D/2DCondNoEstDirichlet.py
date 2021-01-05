from funciones import hdf5
import numpy as np
import sys
import time
from funciones.funciones2D import iterationTime2D, boundary_cond_dirichtlet, GrafContornos, GrafError, Graf3D,Animacion_Contorno,Animacion_Superficie
# =============================================================================
# Especificación de apertura correcta del programa
# =============================================================================
if __name__ == "__main__":
    
    try:
        in_file_name = sys.argv[1]; out_file_name = sys.argv[2]
    except:
        mensaje = """ Error: La ejecucion de este programa requiere de 2 argumentos.
        Ejecucion correcta: python {} entrada salida
        donde "entrada" es el nombre de un archivo que contiene los
        datos del problema :  este se puede generar con el programa hdf5.py.
        El nombre "salida" se usa para almacenar la solucion del problema.

        Por ejemplo: python {} ENTRADA SALIDA""".format(__file__,__file__)


        print(mensaje)
        sys.exit(1)

# =============================================================================
# Lectura de parámetros del archivo de entrada
# =============================================================================
    
    Datos=hdf5.leerParametros(in_file_name,'ax','ay','bx','by','Nx','Ny',\
        'Tx1','Tx2','Ty1','Ty2','ht','Tolerancia','Tmax','ht','Tini','kappa_x','kappa_y','fuente')
    
# =============================================================================
# Imprimimos los parametros y los evaluamos
# =============================================================================
    print('Parametros')
    print('-'*10)
    for key,val in Datos.items(): #ciclo evaluador
        print(key,'=',val)
        print('-'*10)
        exec(key + '=val')      

    hx = (bx-ax)/(Nx+1)  # espaciamiento entre nodos dirección X
    hy = (by-ay)/(Ny+1)  # espaciamiento entre nodos dirección Y

# =============================================================================
# Definimos las condiciones de estabilidad
# =============================================================================

    if kappa_x*ht/hx**2+kappa_y*ht/hy**2<0.5: #ciclo para estabilidad
        pass
    else:

        print('Metodo inestable se cambia el paso en tiempo para cumplir estabilidad')
        print('ht : ',ht,'--->',0.5*(kappa_x/hx**2+kappa_y/hy**2)**-1)
        ht=0.5*(kappa_x/hx**2+kappa_y/hy**2)**-1
    
    Nt=int(Tmax/ht) #Número de pasos en tiempo
    print('Nt = ',Nt)
    print('hx = ',hx)
    print('hy = ',hy)
    
# =============================================================================
# Generación del dominio
# =============================================================================
    x = np.linspace(ax,bx,Nx+2)  #Dominio en dirección X
    y = np.linspace(ay,by,Ny+2)  #Dominio en dirección X
    xg, yg = np.meshgrid(x,y)    # Generación de malla

# =============================================================================
# Definición de la mátriz del sistema
# =============================================================================

    #Aplicamos la condicion inicial a la matriz
    u=np.ones([Ny+2,Nx+2])*Tini
    #Aplicamos las condiciones de frontera tipo dirichlet
    u= boundary_cond_dirichtlet(u,Tx1,Tx2,Ty1,Ty2)

    #Este es el lado derecho de la ecuación, que contiene la condición inicial
    q=np.ones_like(u)*fuente
    # q[5,5]=200500

    errores=[]
    solucion=np.empty([50000,Ny+2,Nx+2],dtype=np.float16)
    Zcambio=[]
    Zcambio.append(u)
    for n in range(Nt+1): #Ciclo para resolver en el espacio
        solucion[n,:,:]=u
        u,error=iterationTime2D(u,q,hx,hy,ht,kappa_x,kappa_y)
        errores.append(error)
        Zcambio.append(u)
        if error < Tolerancia: #Ciclo para Tolerancia
            print('Iteracion terminada con ',n,' pasos')
            break
    
    solucion=solucion[:n+1,:,:]
# =============================================================================
# Preguntar al usuario si desea guardar la solución
# =============================================================================    
    answer=input('Quieres guardar la solucion para generar una animación?  [y/n]\n')
    if answer =='y':
        Datos['solucion_animada']=solucion #Conservar la solución
        Datos['Zcambio']=Zcambio
        # Animacion_Contorno(xg,yg,Zcambio,n)
        # Animacion_Superficie(xg,yg,Zcambio,n)
    elif answer =='n':
        pass
    else:
        pass
# =============================================================================
# Conservar los datos resultantes del calculo en formato hdf5 
# =============================================================================  
    Datos['error']=errores
    Datos['xg']=xg
    Datos['yg']=yg
    Datos['solucion']=u
    Datos['Nt']=n
    Datos['ht']=ht
    
    hdf5.saveParametros(out_file_name,Datos)
    
# =============================================================================
# Generación de la gráfica
# =============================================================================

    #figura 1: mapa de contornos
    GrafContornos(xg, yg, u, 8, 0.75, 'inferno')
    
    #figura 2: superficie
    Graf3D(xg, yg, u, 'inferno')
    
    #figura 3: línea
    GrafError(errores,'C0-')
   
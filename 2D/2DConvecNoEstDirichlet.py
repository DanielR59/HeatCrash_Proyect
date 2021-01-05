from funciones import hdf5
import numpy as np
import sys
from funciones.funciones2D import iterationtimeConvNoEst2D, boundary_cond_dirichtlet, GrafContornos, Graf3D, GrafError,Animacion_Contorno,Animacion_Superficie
# =============================================================================
# Especificación de apertura correcta del programa
# =============================================================================
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

# =============================================================================
# Lectura de parámetros del archivo de entrada
# =============================================================================
    
    
    #leemos los parametros
    Datos=hdf5.leerParametros(in_file_name,'ax','ay','bx','by','Nx','Ny',\
        'Tx1','Tx2','Ty1','Ty2','Tolerancia','Tmax','ht','Tini','kappa_x','kappa_y','fuente','c_p','rho','vel_x','vel_y')
    
# =============================================================================
# Imprimimos los parametros y los evaluamos
# =============================================================================
    print('Parametros')
    print('-'*10)
    for key,val in Datos.items():
        print(key,'=',val)
        print('-'*10)
        exec(key + '=val')
        

    hx = (bx-ax)/(Nx+1)
    hy = (by-ay)/(Ny+1)
    
    

# =============================================================================
# Definimos las condiciones de estabilidad
# =============================================================================
    


    if kappa_x*ht/hx**2+kappa_y*ht/hy**2<0.5:
        pass
    else:

        print('Metodo inestable se cambia el paso en tiempo para cumplir estabilidad')
        print('ht : ',ht,'--->',0.5*(kappa_x/hx**2+kappa_y/hy**2)**-1)
        ht=0.5*(kappa_x/hx**2+kappa_y/hy**2)**-1
        
        


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
    u= boundary_cond_dirichtlet(u,Tx1,Tx2,Ty1,Ty2)


    f=np.ones_like(u)*fuente
    # q[5,5]=200500
    # Definicion de constantes para el calculo
    alfa_x=c_p*rho*vel_x
    alfa_y=c_p*rho*vel_y

    errores=[]
    solucion=np.empty([50000,Ny+2,Nx+2],dtype=np.float16)
    Zcambio=[]
    Zcambio.append(u)
    for n in range(Nt+1):
        solucion[n,:,:]=u
        u,error=iterationtimeConvNoEst2D(u,f,alfa_x,alfa_y,kappa_x,kappa_y,hx,hy,ht)
        errores.append(error)
        Zcambio.append(u)
        if error < Tolerancia:
            print('Iteracion terminada con ',n,' pasos')
            break
    
    solucion=solucion[:n+1,:,:]
    
    answer=input('Quieres guardar la solucion para generar una animación?  [y/n]\n')
    if answer =='y':
        Datos['solucion_animada']=solucion
        Animacion_Contorno(xg,yg,Zcambio,n)
        Animacion_Superficie(xg,yg,Zcambio,n)    
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

# =============================================================================
# Generación de la gráfica
# =============================================================================
    #figura 1: mapa de contornos
    GrafContornos(xg, yg, u, 8, 0.75, 'inferno')
    
    #figura 2: superficie
    Graf3D(xg, yg, u, 'inferno')

    #figura 3: error
    GrafError(errores,'C2-')
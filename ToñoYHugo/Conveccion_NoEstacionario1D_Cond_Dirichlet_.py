#Solución de calor no estacionario conveccion

import numpy as np
import matplotlib.pyplot as plt
import time
import Funciones_No_estacionario as fun

# Programa principal
sel = int(input('''
    CONVECCIÓN 1D - NO ESTACIONARIO

    ¿Cómo deseas ingresar los datos?

    1 - Leer datos desde un archivo
    2 - Ingresar los datos manualmente

    >> '''))


a,b,N,Ta,Tb,K,S,ht,Tmax,v,Tol = fun.lectura_variables(sel,"Datos_Conv_No_EstD.txt") #Ingreso de variables necesarias

h,x,lar,Nt,r,p = fun.constantes(a,b,N,ht,Tmax,K,v) #Cálculo de las demas variables
print("Variables:", "\nAncho de la malla =",h, "\nLongitud de la barra =",lar, "\nPasos de tiempo =",Nt, "\nr =",r, "\np =",p)

# Discretización del dominio
x = fun.dominio(a,b,N)
#Creacion del vector auxiliar
u = fun.vector_aux(Ta,Tb,N,0)

# Este es lado derecho de la ecuación, que contiene la condicion inicial.
# Es decir la solucion en el paso 0. Por eso hacemos una copia de u.
f = np.copy(u[1:N+1])
uold = np.copy(u)

# Construccion de la matriz
A = fun.matriz_Diagonal_Conv_Noest(N,r,p)

# Definición del vector del error
error = []

#Solución Analítica del problema
xa, ua = fun.sol_Analitica(a,b,K,Tmax,N,v)

# Ciclo en el tiempo, desde 0 hasta Nt+1
for n in range(1,Nt+1):
    f[0] += (p+r)*Ta
    f[N-1] += -(p-r)*Tb
    # Solución al sistema lineal
    u[1:N+1] = np.linalg.solve(A,f) 
    
    # Cálculo del error y modificación del vector
    err = np.sqrt(h) * np.linalg.norm(uold-u)
    error.append(err)
    
    if (n % 100==0):
        etiqueta = 'Sol. Num. step = {}'.format(n*ht)
        plt.plot(x, u, '.--', label=etiqueta)
            
    # Actualizacion de la solucion para dar el siguiente paso
    f = np.copy(u[1:N+1])
    uold = np.copy(u)
    
    if (err < Tol):
        break

# Llamado a funciones para la graficación
fun.graficas(xa,ua,'Solución Analítica y Numérica')
fun.grafica_error(error)

# Llamado a función para guardar los datos
fun.escritura(xa,ua)

# Solucion de calor no estacionario conveccion

import Funciones_No_estacionario as fun
import numpy as np
import matplotlib.pyplot as plt
import time

# Programa principal
sel = int(input('''
    CONVECCIÓN 1D - NO ESTACIONARIO

    ¿Cómo deseas ingresar los datos?

    1 - Leer datos desde un archivo
    2 - Ingresar los datos manualmente

    >> '''))

a,b,N,Ta,Tb,K,S,ht,Tmax,v,Tol = fun.lectura_variables(sel,"Datos_Conv_No_EstD.txt")

# Cálculo de Constantes
h,x,lar,Nt,r,p = fun.constantes(a,b,N,ht,Tmax,K,v)
print("\tVariables:", "\n\tAncho de la malla =",h, "\n\tLongitud de la barra =",lar, "\n\tPasos de tiempo =",Nt, "\n\tr =",r, "\n\tp =",p)

x = np.linspace(a,b,N+2)

# Vector que contiene las fuentes
q = np.ones(N+2) * (-S)

#Creacion del vector auxiliar
u = fun.vector_aux(Ta,Tb,N,q)

#Lado derecho de la ecuación:
f = np.copy(u[1:N+1])
uold = np.copy(u)

# Construccion de la matriz
A = fun.matriz_Diagonal_Conv_Noest(N,r,p)

#Graficación de la condición inicial
plt.plot(x,u,'--k',label='Inicial')

tolerancia=Tol
error = []

#Solución Analítica del problema
xa, ua = fun.sol_Analitica(a,b,K,Tmax,N,v)

# Ciclo en el tiempo, desde 0 hasta Nt+1
for n in range(1,Nt+1):
    
    f[0] += (p+r)*Ta
    f[N-1] += -(p-r)*Tb
    u[1:N+1] = np.linalg.solve(A,f) # Sol. del sistema lineal
    
    err = np.sqrt(h) * np.linalg.norm(uold-u)
    error.append(err)
    
    #print("n = ", n, ' Error = %12.10g' % err)
    if (n % 100==0):
        etiqueta = 'Sol. Numer. Paso = {}'.format(n*ht)
        plt.plot(x, u, '.--', label=etiqueta)
        
        
    # Actualizacion de la solucion para dar el siguiente paso
    f = np.copy(u[1:N+1])
    uold = np.copy(u)
    

    if (err < tolerancia):
        break


fun.graficas(xa,ua,'Solución Analítica y Numérica')
fun.grafica_error(error)

import hdf5
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation

def plotheatmap(xg,yg,u_k, Nt,ht,minimo,maximo):
    # Clear the current plot figure
        plt.clf()
        plt.title(f"Temperatura en t = {Nt*ht:.3f} s")
        plt.xlabel("x")
        plt.ylabel("y")
        
    # This is to plot u_k (u at time-step k)
        plt.pcolormesh(xg,yg,u_k, cmap=plt.cm.jet, vmin=minimo, vmax=maximo,shading='nearest')
        plt.colorbar()
    
        return plt

if __name__ == "__main__":
    
    try:
        in_file_name = sys.argv[1]; 
    except:
        mensaje = """ Error: La ejecucion de este programa requiere de 2 argumentos.
        Ejecucion correcta: python 2D_Poisson_01.py entrada salida
        donde "entrada" es el nombre de un archivo que contiene los
        datos del problema :  este se puede generar con el programa hdf5.py.
        El nombre "salida" se usa para almacenar la solucion_animada del problema.

        Por ejemplo: python 2D_Poisson_01.py SALIDA""" 

        print(mensaje)
        sys.exit(1)

    
    #leemos los parametros
    Datos=hdf5.leerParametros(in_file_name,'Nt','ht','xg','yg','solucion_animada')
    for key,val in Datos.items():

        exec(key + '=val')
    
    if Nt>=20000:
        step=4000
    elif 10000<=Nt<20000:
        step = 2000
    elif 6000<Nt<10000:
        step=200
    else:
        step=50
    
    def animate(Nt):
        minimo=np.amin(solucion_animada)
        maximo=np.amax(solucion_animada)
        plotheatmap(xg,yg,solucion_animada[Nt], Nt,ht,minimo,maximo)
    print("Empezando a generar animacion")
    print("Este proceso puede tardar varios minutos dependiendo del numero de soluciones Nt \n ve y echate un refresquito")
    fig = plt.figure()
    anim = FuncAnimation(fig, animate,frames=range(0,Nt+1,step), interval=500, repeat=False)
    anim.save("heat_equation_solution.gif")
    print("Animacion lista :D \nBuscala en tu carpeta")
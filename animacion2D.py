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
        in_file_name = sys.argv[1]; out_file_name = sys.argv[2]
    except:
        mensaje = """ Error: La ejecucion de este programa requiere de 2 argumentos.
        Ejecucion correcta: python {} entrada salida
        donde "entrada" es el nombre de un archivo que contiene los
        datos resultados de un problema 2D que involucre el tiempo.
        El nombre "salida" es el nombre y tipo de archivo que se generara de la animación.
        Este puede ser tanto .gif como .mp4
        NOTA: Si no tienes ffmpeg no podras guardar como .mp4

        Por ejemplo: python {} ENTRADA SALIDA.gif""".format(__file__,__file__) 

        print(mensaje)
        sys.exit(1)

    
    #leemos los parametros
    Datos=hdf5.leerParametros(in_file_name,'Nt','ht','xg','yg','solucion_animada')
    for key,val in Datos.items():

        exec(key + '=val')
    #Sleccionamos un paso adecuado
    if Nt>=20000:
        step = int(0.05/ht)
    elif 10000<=Nt<20000:
        step = int(0.02/ht)
    elif 6000<Nt<10000:
        step = int(0.005/ht)
    else:
        step = 50
    
    #Checamos el nombre del archivo de salida
    if out_file_name.endswith('.gif'):
        pass
    elif out_file_name.endswith('.mp4'):
        pass
    else:
        out_file_name+='.gif'


    #Parametros de la animacion
    minimo=np.amin(solucion_animada)
    maximo=np.amax(solucion_animada)
    def animate(Nt):        
        plotheatmap(xg,yg,solucion_animada[Nt], Nt,ht,minimo,maximo)

    #Empezamos a generar la animacion    
    print("Empezando a generar animacion")
    print("Este proceso puede tardar varios minutos dependiendo del numero de soluciones Nt \nVe y echate un refresquito")
    fig = plt.figure()
    anim = FuncAnimation(fig, animate,frames=range(0,Nt+1,step), interval=500, repeat=False)
    anim.save(out_file_name)
    
    print("Animacion lista :D \nBuscala en tu carpeta")
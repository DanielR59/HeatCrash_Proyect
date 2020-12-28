import h5py
import numpy as np
import sys
# =============================================================================
# Función para la lectura de parámetros usando la herramienta diccionarios
# =============================================================================
def leerParametros(archivoHDF5,*parametros,diccionario={}):
    with h5py.File(archivoHDF5,'r') as f:
        for parametro in parametros:
            valor=f.get(parametro)
            if valor==None:
                print('Parametro :', parametro,' no encontrado')
                continue

            if valor.dtype == int:
                diccionario[parametro]=int(np.array(valor))

            elif valor.dtype == float:
                try:
                    diccionario[parametro]=float(np.array(valor))
                except TypeError:
                    diccionario[parametro]=np.array(valor)
            elif valor.dtype == np.float16:
                diccionario[parametro]=np.array(valor,dtype=np.float16)
            elif valor.dtype == object:
                diccionario[parametro]=np.array2string(valor)
        return diccionario

# =============================================================================
# Función para guardar los parámetros en el archivo HDF5
# =============================================================================

def saveParametros(archivoHDF5,diccionario):


    with h5py.File(archivoHDF5,'w') as f:
        for key,values in diccionario.items():
            f[key]=values



if __name__ == "__main__":
    try:
        in_file_name=sys.argv[1]
    except:
        mensaje = """ Error: La ejecucion de este programa requiere de 1 argumento.
        Ejecucion correcta: python hdf5.py ENTRADA
        donde el nombre "ENTRADA" es el nombre del archivo donde
        se almacenan los datos del problema.
        Los datos que se guardan en el archivo ENTRADA son los que se encuentran en
        este archivo hdf5.py

        Por ejemplo: python hdf5.py ENTRADA"""

        print(mensaje)
        sys.exit(1)

    #Generamos parametros
    archivo=in_file_name
    Datos={

        'ax' : 0, #Punto inicial en dirección x
        'bx' : 1,  #Punto final en dirección x
        'ay' : 0,  #Punto inicial en dirección y
        'by' : 1, #Punto final en dirección y
        'Nx' : 60, # Número de nodos en x
        'Ny' : 60, # Número de nodos en y
        'Tx1' : 40, # Condición de frontera 1 en x
        'Tx2' : -20, #Condición de frontera 2 en x
        'Ty1' : 0, # Condición de frontera 1 en y
        'Ty2' : -300, # Condición de frontera 2 en x
        'Tini' : 0, # Condición inicial
        'ht' : 0.01, # Espaciamiento en tiempo
        'Tmax' : 1, # Tiempo máximo
        'Tolerancia' : 1E-4, # Tolerancia

    }

    saveParametros(archivo,Datos)

import h5py
import numpy as np 
import sys
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
            elif valor.dtype == object:
                diccionario[parametro]=np.array2string(valor) 
        return diccionario


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
        'ax' : 0,
        'bx' : 1,
        'ay' : 0,
        'by' : 1,
        'Nx' : 50,
        'Ny' : 50,
        'Tx1' : 100,
        'Tx2' : 50,
        'Ty1' : 0,
        'Ty2' : 200,
        'Tini' : 10,
        'ht' : 0.01,
        'Tmax' : 1,
        'Tolerancia' : 1E-4,

    }

    saveParametros(archivo,Datos)


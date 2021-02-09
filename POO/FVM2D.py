

from Mesh2D import Mesh2D
from Coefficients2D import Coefficients2D
from Diffusion2D import Diffusion2D
from Advection2D import Advection2D
from Temporal2D import Temporal2D
from Matrix2D import Matrix2D


def crono(f):
    """
    Regresa el tiempo que toma en ejecutarse la funcion.
    """
    def eTime(A,b):
        t1 = time.time()
        f(A,b)
        t2 = time.time()
        return 'Elapsed time: ' + str((t2 - t1)) + "\n"
    return eTime

def decorate(f):
    def nicePrint(**kargs):
        line = '-' * 70
        print('.'+ line + '.')
        print('|{:^70}|'.format('NoNacos : Numerical Objects for Natural Convection Systems'))
        print('.'+ line + '.')
        print('|{:^70}|'.format(' Ver. 0.1, Author LMCS, 2018, [GNU GPL License V3]'))
        print('.'+ line + '.')
        f(**kargs)
        print('.'+ line + '.')
    return nicePrint

@decorate
def printData(**kargs):
    for (key,value) in kargs.items():
        if (type(value) == str):
            print('|{:^70}|'.format('{0:>15s} = {1:11s}'.format(key, value)))
        elif (type(value) == int):
            print('|{:^70}|'.format('{0:>15s} = {1:<11d}'.format(key, value)))            
        else:
            print('|{:^70}|'.format('{0:>15s} = {1:10.5e}'.format(key, value)))

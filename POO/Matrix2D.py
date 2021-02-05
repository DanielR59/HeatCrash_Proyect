



import numpy as np

class Matrix2D():

    def __init__(self, nvx = None, nvy = None):
        

        self.Nx = nvx -2
        self.Ny = nvy -2
        self.N = (nvx-2)*(nvy-2)
        self.A = np.eye(self.N)

    def build(self, coefficients = None):

        aP = coefficients.aP
        aE = coefficients.aE
        aW = coefficients.aW
        aN = coefficients.aN
        aS = coefficients.aS
        Nx = self.Nx
        Ny = self.Ny
        A = self.A
        multiples_auxiliares=np.array([i for i in range(Nx,Nx*Ny,Nx)]) #No preguntes solo gozalo
        multiples_auxiliares-=1 #Porque python empieza en 0 

        aP_aux = aP[1:-1,1:-1].ravel()
        aE_aux = aE[1:-1,1:-1].ravel()
        aW_aux = aW[1:-1,1:-1].ravel()
        aS_aux = aS[1:-1,1:-1].ravel()
        aN_aux = aN[1:-1,1:-1].ravel()
        
        for i in range(0,Nx*Ny-1):
            A[i,i] = aP_aux[i]
            if (i not in multiples_auxiliares):
                A[i+1,i] = aW_aux[i]
                A[i,i+1] = aE_aux[i]

        for i in range(0,Nx*Ny-Nx): #Se llenan los valores de la matriz 
            A[Nx+i,i]=aN_aux[i] 
            A[i,Nx+i] = aS_aux[i]



        


if __name__ == '__main__':

    Matriz = Matrix2D(5,5)
    
    from Diffusion2D import Diffusion2D

    df1 = Diffusion2D(5,5,1,1,0.25)
    df1.alloc()
    df1.calcCoef()
    df1.bcDirichlet('LEFT_WALL', 2)
    df1.bcDirichlet('RIGHT_WALL', 1)
    
    print(df1.aP)

    Matriz.build(df1)
    print('-'*10)
    print(Matriz.A)


        
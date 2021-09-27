
import numpy as np

from Coefficients2D import Coefficients2D

class Temporal2D(Coefficients2D):

    def __init__(self,nvx : int = None,nvy : int = None, rho : float = None, dx : float = None,dy : float = None, dt : float = None):
        super().__init__(nvx,nvy)
        self.nvx = nvx
        self.nvy = nvy
        self.rho = rho
        self.dx = dx
        self.dy = dy
        self.dt = dt

    def calcCoef(self, phi_old):
        aP = self.aP
        Su = self.Su
        rho = self.rho
        dxy_dt = self.dx*self.dy/self.dt

        # Sustitución del for loop
        aP[1:-1,1:-1] += rho * (dxy_dt)
        Su[1:-1,1:-1] += phi_old[1:-1,1:-1] * (dxy_dt)

        # for i in range(1,self.nvx-1):
        #     for j in range(1, self.nvy-1):
        #         #DUDA
        #         aP[i,j] += rho * (dxy_dt)
        #         Su[i,j] += phi_old[i,j]*(dxy_dt)

if __name__ == '__main__':

    nvx = 6
    nvy = 6
    
    aux = np.random.rand(nvx,nvy)

    phi_old = np.sin(aux)

    Tf = Temporal2D(6,6,1,1,1,1)
    Tf.alloc()
    Tf.calcCoef(phi_old)
    print(Tf.Su)






    
    
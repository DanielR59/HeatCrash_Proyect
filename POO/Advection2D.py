import numpy as np
from Coefficients2D import Coefficients2D 


class Advection2D(Coefficients2D):

    def __init__(self, nvx=None, nvy=None, deltaX=None, deltaY=None, rho = None):
        super().__init__(nvx=nvx, nvy=nvy, deltaX=deltaX, deltaY=deltaY)
    
        self.nvx=nvx
        self.nvy = nvy
        self.rho = rho
        self.dx = deltaX
        self.dy = deltaY
        self.u = np.zeros([nvx-1,nvy-1])
        

    def setU(self, u):
        self.u+=u
        

    def calcCoef(self):

        aE = self.aE
        aW = self.aW
        aP = self.aP
        aN = self.aN
        aS = self.aS
        u = self.u
        rho = self.rho
        nvx = self.nvx
        nvy = self.nvy

    

        CE = - rho * u[1:,:]*0.5
        CW = rho * u[:-1,:]*0.5
        CN = -rho * u[:,:-1]*0.5
        CS = rho * u[:,1:]*0.5

        for i in range(1,nvx-1):
            for j in range(1,nvy-1):
                CE = - rho * u[i,j]*0.5
                CW = rho * u[i-1,j]*0.5
                CN = -rho * u[i,j]*0.5
                CS = rho * u[i,j-1]*0.5
                aE[i,j]+=CE
                aW[i,j]+=CW
                aS[i,j]+=CS
                aN[i,j]+=CN


                aP[i,j] += CE+CN+CW+CS+rho*(u[i,j]-u[i-1,j])+rho* (u[i,j]-u[i,j-1])



        print(aP)





if __name__== '__main__':

    Holi = Advection2D(4,4,1,1,2)
    Holi.alloc()
    Holi.setU(np.eye(3)*5.9)
    Holi.calcCoef()

    Holi.bcDirichlet('TOP_WALL',5)



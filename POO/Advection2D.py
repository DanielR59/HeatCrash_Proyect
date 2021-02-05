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
        self.ux = np.zeros([nvx-1,nvy-1])
        self.uy = np.zeros([nvx-1,nvy-1])


    def setUx(self, u):
        self.ux+=u
        
    def setUy(self, u):
        self.uy+=u

    def calcCoef(self):

        aE = self.aE
        aW = self.aW
        aP = self.aP
        aN = self.aN
        aS = self.aS
        ux = self.ux
        uy = self.uy
        rho = self.rho
        nvx = self.nvx
        nvy = self.nvy

        
        for i in range(1,nvx-1):
            for j in range(1,nvy-1):
                CE = - rho * ux[i,j]*0.5
                CW = rho * ux[i-1,j]*0.5
                CN = -rho * uy[i,j]*0.5
                CS = rho * uy[i,j-1]*0.5
                aE[i,j]+=CE
                aW[i,j]+=CW
                aS[i,j]+=CS
                aN[i,j]+=CN


                aP[i,j] += CE+CN+CW+CS+rho*(ux[i,j]-ux[i-1,j])+rho* (uy[i,j]-uy[i,j-1])








if __name__== '__main__':

    Holi = Advection2D(4,4,1,1,2)
    Holi.alloc()
    Holi.setU(np.eye(3)*5.9)
    Holi.calcCoef()

    Holi.bcDirichlet('TOP_WALL',5)



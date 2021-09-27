
from Coefficients2D import Coefficients2D 
import numpy as np

class Diffusion2D(Coefficients2D):


    def __init__(self, nvx : int =None, nvy:int=None, deltaX:float=None, deltaY:float=None, Gamma:float = None):
        super().__init__(nvx=nvx, nvy=nvy, deltaX=deltaX, deltaY=deltaY)
    
        self.nvx=nvx
        self.nvy = nvy
        self.Gamma = Gamma
        self.dx = deltaX
        self.dy = deltaY
        

    def calcCoef(self):

        aE = self.aE
        aW = self.aW
        aP = self.aP
        aN = self.aN
        aS = self.aS
        Gamma =self.Gamma 
        deltaX=self.dx
        deltaY = self.dy

        aE+=Gamma*deltaY/deltaX
        aW+=Gamma*deltaY/deltaX
        aN+=Gamma*deltaX/deltaY
        aS+=Gamma*deltaX/deltaY

        aP += aE+aN+aS+aW 

if __name__ == '__main__':

    df = Diffusion2D(5,5,1,1,0.6)
    df.alloc()
    df.calcCoef()
    df.bcDirichlet('TOP_WALL',5)
    print(df.aP)

import numpy as np
#CHECAR DUDAS
class Coefficients2D():



    aP = None
    aE = None
    aW = None
    aN = None
    aS = None
    Su = None
    nvx = None
    nvy = None
    deltaX = None
    deltaY = None

    def __init__(self, nvx = None, nvy = None, deltaX = None, deltaY = None):
        Coefficients2D.nvx=nvx
        Coefficients2D.nvy=nvy
        Coefficients2D.deltaX=deltaX
        Coefficients2D.deltaY=deltaY
    
    @staticmethod
    def alloc():
        nvx = Coefficients2D.nvx
        nvy = Coefficients2D.nvy

        Coefficients2D.aP = np.zeros((nvx,nvy))
        Coefficients2D.aE = np.zeros((nvx,nvy))
        Coefficients2D.aW = np.zeros((nvx,nvy))
        Coefficients2D.aS = np.zeros((nvx,nvy))
        Coefficients2D.aN = np.zeros((nvx,nvy))
        Coefficients2D.Su = np.zeros((nvx,nvy))
    @staticmethod
    #CHECAR TEORIA
    def bcDirichlet(wall, phi):
        aP = Coefficients2D.aP
        aE = Coefficients2D.aE
        aW = Coefficients2D.aW
        Su = Coefficients2D.Su
        aN = Coefficients2D.aN
        aS = Coefficients2D.aS

        if wall == 'LEFT_WALL':
            aP[:,1] += aW[:,1]
            Su[:,1] += 2 * aW[:,1] * phi
        elif wall == 'RIGHT_WALL':
            aP[:,-2] += aE[:,-2]
            Su[:,-2] += 2 * aE[:,-2] * phi
        
        elif wall == 'TOP_WALL':
            aP[1,:] += aN[1,:]
            Su[1,:] += 2 * aN[1,:] * phi
        elif wall == 'BOTTOM_WALL':
            aP[-2,:] += aS[-2,:]
            Su[-2,:] += 2 * aS[-2,:] * phi
    @staticmethod
    def bcNeumman(wall, flux):
        aP = Coefficients2D.aP
        aE = Coefficients2D.aE
        aW = Coefficients2D.aW
        aN = Coefficients2D.aN
        aS = Coefficients2D.aS
        Su = Coefficients2D.Su
        dx = Coefficients2D.deltaX
        dy = Coefficients2D.deltaY

        if wall == 'LEFT_WALL':
            aP[:,0] -= aW[:,0]
            Su[:,1] -= aW[:,0] * flux * dx
        elif wall == 'RIGHT_WALL':
            aP[:,-1] -= aE[:,-1]
            Su[:,-2] += aE[:,-1] * flux * dx

        elif wall == 'TOP_WALL':
            aP[0,:] -= aN[0,:]
            Su[1,:] -= aN[0,:] * flux * dy
            
        elif wall == 'BOTTOM_WALL':
            aP[-1,:] -= aE[-1,:]
            Su[-2,:] += aE[-1,:] * flux * dy
#DUDA
    def setSu(self, q):
        Su = Coefficients2D.Su
        dx = Coefficients2D.deltaX
        Su += q * dx
    #Aqui igual   
    def cleanCoefficients(self):
        Coefficients2D.aP[:,:] = 0
        Coefficients2D.aE[:,:] = 0
        Coefficients2D.aW[:,:] = 0
        Coefficients2D.aS[:,:] = 0
        Coefficients2D.aN[:,:] = 0
        Coefficients2D.Su[:,:] = 0
        

if __name__ == '__main__':


    Algo=Coefficients2D(8,8, 0.1, 0.1)

    Algo.alloc()
    Algo.aP[:,:] = 5
    Algo.aW[:,:] = 2
    Algo.aE[:,:] = 7
    Algo.aN[:,:] = 4
    Algo.aS[:,:] = 3
    Algo.setSu(30)
    print(Algo.aP)
    Algo.bcDirichlet('LEFT_WALL',20)
    Algo.bcDirichlet('RIGHT_WALL',20)
    Algo.bcDirichlet('TOP_WALL',20)
    Algo.bcDirichlet('BOTTOM_WALL',20)


    print(Algo.aP)
    print(Algo.Su)
        


import numpy as np


class Mesh2D():

    def __init__(self, nodesX = None, nodesY = None, volumesX = None, volumesY = None, lengthX = None, lengthY = None):

        self.nodesX = nodesX
        self.nodesY = nodesY
        self.volumesX = volumesX
        self.volumesY = volumesY
        self.lengthX = lengthX
        self.lengthY = lengthY
        self.calcDeltaX()
        self.calcDeltaY()


    def calcDeltaX(self):
            self.deltaX = self.lengthX / (self.nodesX - 1)
            
    def calcDeltaY(self):
            self.deltaY = self.lengthY / (self.nodesY - 1)
    def createMesh(self):
        first_volumeX = self.deltaX / 2
        final_volumeX = self.lengthX - first_volumeX
        first_volumeY = self.deltaY / 2
        final_volumeY = self.lengthY - first_volumeY

        self.__x = np.zeros(self.volumesX)
        self.__y = np.zeros(self.volumesY)
        
        self.__x[1:-1] = np.linspace(first_volumeX,final_volumeX,self.volumesX-2)
    
        self.__y[1:-1] = np.linspace(first_volumeY,final_volumeY,self.volumesY-2)

        self.__x[-1] =self.lengthX
        self.__y[-1] =self.lengthY

        self.X,self.Y = np.meshgrid(self.__x,self.__y)
        print('Malla X\n',self.X,'Malla Y\n',self.Y)
        return self.X,self.Y

if __name__ == '__main__':
    objeto=Mesh2D(3,3,4,3,1.5,6.2)
    print(objeto.createMesh())


    
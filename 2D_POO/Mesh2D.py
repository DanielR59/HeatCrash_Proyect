import numpy as np


class Mesh2D:

    def __init__(self, nodesX : int = None, nodesY : int = None, volumesX : int = None, volumesY : int = None, lengthX : float = None, lengthY : float = None):

        self.nodesX = nodesX
        self.nodesY = nodesY
        self.volumesX = volumesX
        self.volumesY = volumesY
        self.lengthX = lengthX
        self.lengthY = lengthY
        self.calcDeltaX()
        self.calcDeltaY()
        self.adjustNodesVolumesX(nodesX, volumesX)
        self.adjustNodesVolumesY(nodesY, volumesY)


    def calcDeltaX(self):
        self.dx = self.lengthX / (self.nodesX - 1)
            
    def calcDeltaY(self):
        self.dy = self.lengthY / (self.nodesY - 1)

    def adjustNodesVolumesX(self,nodes : int,volumes : int):
        if nodes:
            self.volumesX = self.nodesX + 1
        if volumes:
            self.nodesX = self.volumesX - 1        
    
    def adjustNodesVolumesY(self,nodes : int,volumes : int):
        if nodes:
            self.volumesY = self.nodesY + 1
        if volumes:
            self.nodesY = self.volumesY - 1        
        


    def createMesh(self):
        first_volumeX = self.dx / 2
        final_volumeX = self.lengthX - first_volumeX
        first_volumeY = self.dy / 2
        final_volumeY = self.lengthY - first_volumeY

        self.__x = np.zeros(self.volumesX)
        self.__y = np.zeros(self.volumesY)
        
        self.__x[1:-1] = np.linspace(first_volumeX,final_volumeX,self.volumesX-2)
    
        self.__y[1:-1] = np.linspace(first_volumeY,final_volumeY,self.volumesY-2)

        self.__x[-1] =self.lengthX
        self.__y[-1] =self.lengthY

        self.X,self.Y = np.meshgrid(self.__x,self.__y)

if __name__ == '__main__':
    objeto=Mesh2D(3,3,4,3,1.5,6.2)
    objeto.createMesh()
    print(objeto.X)

    
    
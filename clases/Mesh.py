# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 14:33:25 2021

@author: josuelg

"""
import numpy as np

class Mesh():
    
    def __init__(self, Lx, Nx, Ly=None, Ny= None):
        self.Lx=Lx
        self.Nx=Nx
        self.Ly=Ly
        self.Ny=Ny
        
        print('Lx = {}, NX = {}'. format(self.Lx,self.Nx))
        if (self.Ly!=None and self.Ny!=None):
            print('Ly = {}, Ny = {}'. format(self.Ly,self.Ny))

    def Delta(self):
        self.hx=self.Lx/(self.Nx-1)
        print('hx = {}'. format(self.hx))
        if (self.Ly!=None and self.Ny!=None):
            self.hy=self.Ly/(self.Ny-1)
            print('hy = {}'. format(self.hy))
            return self.hx, self.hy
        return self.hx
    
    def createMesh(self):
        if(self.Ly!=None and self.Ny!=None):
            self.malla=np.zeros((self.Nx,self.Ny))
            self.x=np.linspace(0,self.Lx,self.Nx)
            self.y=np.linspace(0,self.Ly,self.Ny)
            return self.x, self.y, self.malla
        else:
            self.malla=np.zeros(self.Nx)
            self.x=np.linspace(0,self.Lx,self.Nx)
            return self.x, self.malla

if __name__ == '__main__':
    malla=Mesh(1.0,11,2.0,21)
    hx=malla.Delta()
    x,y,T=malla.createMesh()


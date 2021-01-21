#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 17:32:51 2021

@author: josuelg
"""
import numpy as np

class Matrix():
    
    def __init__(self,n,m):
        self.n=n
        self.m=m
        
    def crea_Matriz(self):
        self.matriz=np.zeros((self.n,self.m))
        
    def A_Dir(self):
        for i in range(self.n):
            for j in range(self.m):
                if (i==j):
                    self.matriz[i][j]=-2
                if(np.fabs(i-j)==1):
                    self.matriz[i][j]=1
        return self.matriz
    
    def A_New(self):
        for i in range(self.n):
            for j in range(self.m):
                if (i==j):
                    self.matriz[i][j]=-2
                if(np.fabs(i-j)==1):
                    self.matriz[i][j]=1
        self.matriz[self.n-1][self.m-2]=-1
        self.matriz[self.n-1][self.m-1]=1
        return self.matriz        

if __name__ == '__main__':
    matriz=Matrix(5,5)
    matriz.crea_Matriz()
    A=matriz.A_New()
    print(A)
    B=matriz.A_Dir()
    print(B)
    
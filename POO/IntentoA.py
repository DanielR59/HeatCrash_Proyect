import numpy as np
Nx=10
Ny = 5
A = np.eye(Nx*Ny)

aP = np.zeros([Nx,Ny])

for i in range(0,Nx*Ny):
    A[i,i]=aP.ravel()[i]


import FVM2D as fvm
import numpy as np
import matplotlib.pyplot as plt



longitudx = 0.5 # meters
longitudy = 0.5
TA = 100 # °C 
TB = 500 # °C 
TC = 100
TD = 500
k  = 1000 # W/m.K
Nx  = 6 # Número de nodos
Ny = 6


malla = fvm.Mesh2D(6,6,lengthX=longitudx, lengthY=longitudy)
malla.calcDeltaX()
print(malla.nodesX,malla.volumesX)


df2 = fvm.Diffusion2D(nvx = malla.volumesX,nvy = malla.volumesY, deltaX= malla.dx, deltaY = malla.dy, Gamma=k)
print(malla.dx,malla.dy)
df2.alloc()
df2.calcCoef()



T = np.zeros([df2.nvx,df2.nvy])
T[:,0] = TA
T_aux=T[1:-1,1:-1].ravel()
Nx=df2.nvx-2
Ny = df2.nvy-2
df2.bcDirichlet('LEFT_WALL',TA)

print(T.size)

A = fvm.Matrix2D(malla.volumesX,malla.volumesY)

A.build(df2)
print(T[1:-1,1:-1].size,A.A.shape,df2.Su[1:-1,1:-1].size)
T_aux = np.linalg.solve(A.A,df2.Su[1:-1,1:-1].ravel())
T_aux.shape = (Nx,Ny)

T[1:-1,1:-1] =T_aux

print(df2.Su)


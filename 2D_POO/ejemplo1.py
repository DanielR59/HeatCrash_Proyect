
import FVM2D as fvm
import numpy as np
import matplotlib.pyplot as plt



longitudx = 1 # meters
longitudy = 1
TA = 100 # °C 
TB = 500 # °C 
TC = 100
TD = 500
k  = 1000 # W/m.K

#Generamos la malla
malla = fvm.Mesh2D(40,40,lengthX=longitudx, lengthY=longitudy)
print(malla.nodesX,malla.volumesX)
malla.createMesh()
print(malla.dx)
#Generamos, almacenamos y calculamos los coeficientes del problema de difusion
df2 = fvm.Diffusion2D(nvx = malla.volumesX,nvy = malla.volumesY, deltaX= malla.dx, deltaY = malla.dy, Gamma=k)
df2.alloc()
df2.calcCoef()
# print(df2.aP,df2.aE,sep='\n')

# df2.setSu(500)

#Generamos el vector que almacenara la solucion
T = np.zeros([df2.nvx,df2.nvy])
#Le aplicamos las condiciones de frontera a los coeficientes y a el vector solucion

T[:,0] = TA
T[:,-1] = TD
T[0,:] = TC
T[-1,:] = TB
T_aux=T[1:-1,1:-1].ravel()
Nx=df2.nvx-2
Ny = df2.nvy-2
df2.bcDirichlet('LEFT_WALL',TA)
df2.bcNeumman('TOP_WALL',TC)
df2.bcDirichlet('RIGHT_WALL',TD)
df2.bcNeumman('BOTTOM_WALL',TB)
#
#GENERAMOS Y CALCULAMOS LA MATRIZ A PARTIR DE LOS COEFICIENTES
A = fvm.Matrix2D(malla.volumesX,malla.volumesY)
A.build(df2)
print(A.A)
# print(np.linalg.det(A.A))
# print(np.linalg.inv(A.A)) #Probablemente lo que está mal es A

print(T[1:-1,1:-1].size,A.A.shape,df2.Su[1:-1,1:-1].size)
print(df2.Su[1:-1,1:-1].flatten())

#Solucionamos el sistema lineal
T_aux = np.linalg.solve(A.A,df2.Su[1:-1,1:-1].flatten())
T_aux.shape = (Nx,Ny)

T[1:-1,1:-1] =T_aux
# print(T)

f1 = plt.figure()
c = plt.contourf(malla.X,malla.Y,T,8, alpha=.75,cmap='inferno')
f1.colorbar(c, shrink=1.0)
plt.show()





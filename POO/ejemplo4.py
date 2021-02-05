
import FVM2D as fvm
import numpy as np
import matplotlib.pyplot as plt

longitudx = 1 # meters
longitudy = 1
TA = 100 # °C 
TB = 500 # °C 
TC = 100
TD = 500
q = 1e+06
rho = 1.0 # kg/m^3
k  = .1 # W/m.K
ux = 2.1 # m/s
uy = 2.1 # m/s

#Generamos la malla
malla = fvm.Mesh2D(40,40,lengthX=longitudx, lengthY=longitudy)
print(malla.nodesX,malla.volumesX)
malla.createMesh()
print(malla.dx)

coef = fvm.Coefficients2D(malla.volumesX,malla.volumesY,malla.dx,malla.dy)
coef.alloc()

dif = fvm.Diffusion2D(malla.volumesX,malla.volumesY,malla.dx,malla.dy,k)

dif.calcCoef()

#  Calculamos los coeficientes de FVM de la Advección

adv = fvm.Advection2D(malla.volumesX,malla.volumesY,malla.dx,malla.dy,rho)

adv.setUx(ux)
adv.setUy(uy)
adv.calcCoef()


#Generamos el vector que almacenara la solucion
T = np.zeros([coef.nvx,coef.nvy])
#Le aplicamos las condiciones de frontera a los coeficientes y a el vector solucion

T[:,0] = TA
T[:,-1] = TD
T[0,:] = TC
T[-1,:] = TB
T_aux=T[1:-1,1:-1].ravel()
Nx=coef.nvx-2
Ny = coef.nvy-2
coef.bcDirichlet('LEFT_WALL',TA)
coef.bcDirichlet('TOP_WALL',TC)
coef.bcDirichlet('RIGHT_WALL',TD)
coef.bcDirichlet('BOTTOM_WALL',TB)
# Se construye el sistema lineal de ecuaciones a partir de los coef. de FVM
#
A = fvm.Matrix2D(malla.volumesX,malla.volumesY)
A.build(coef)
# print(A.A)

#Solucionamos el sistema lineal
T_aux = np.linalg.solve(A.A,coef.Su[1:-1,1:-1].flatten())
T_aux.shape = (Nx,Ny)

T[1:-1,1:-1] =T_aux
# print(T)


f1 = plt.figure()
c = plt.contourf(malla.X,malla.Y,T,8, alpha=.75,cmap='inferno')
f1.colorbar(c, shrink=1.0)
plt.show()





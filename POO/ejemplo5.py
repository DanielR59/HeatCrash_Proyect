

import FVM2D as fvm
import numpy as np
import matplotlib.pyplot as plt

longitudx = 10 # meters
longitudy = 10
TA = 100 # °C 
TB = 200 # °C 
TC = 100
TD = 500
q = 1e+06
rho = 1.0 # kg/m^3
k  = 1000 # W/m.K
ux = .1 # m/s
uy = .1 # m/s

delta_t = 0.002 # Paso de tiempo
steps = 500

malla = fvm.Mesh2D(15,15,lengthX=longitudx, lengthY=longitudy)
print(malla.nodesX,malla.volumesX)
malla.createMesh()
print(malla.dx)

coef = fvm.Coefficients2D(malla.volumesX,malla.volumesY,malla.dx,malla.dy)
coef.alloc()

dif = fvm.Diffusion2D(malla.volumesX,malla.volumesY,malla.dx,malla.dy,k)
adv = fvm.Advection2D(malla.volumesX,malla.volumesY,malla.dx,malla.dy,rho)
tem = fvm.Temporal2D(malla.volumesX,malla.volumesY,rho,malla.dx,malla.dy,delta_t)

T = np.ones([coef.nvx,coef.nvy])*0 #Condicion inicial
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



for i in range(1, steps+1):

    coef.cleanCoefficients()
    dif.calcCoef()
    adv.setUx(ux)
    adv.setUy(uy)
    adv.calcCoef()
    tem.calcCoef(T)

    coef.bcDirichlet('LEFT_WALL',TA)
    coef.bcDirichlet('TOP_WALL',TC)
    coef.bcDirichlet('RIGHT_WALL',TD)
    coef.bcDirichlet('BOTTOM_WALL',TB)

    A = fvm.Matrix2D(malla.volumesX,malla.volumesY)
    A.build(coef)

    T_aux = np.linalg.solve(A.A,coef.Su[1:-1,1:-1].flatten())
    T_aux.shape = (Nx,Ny)
    T[1:-1,1:-1] =T_aux

    


f1 = plt.figure()
c = plt.contourf(malla.X,malla.Y,T,8, alpha=.75,cmap='inferno')
f1.colorbar(c, shrink=1.0)


surf=plt.figure(figsize=(5,4)) 
ax=surf.gca(projection='3d')
s=ax.plot_surface(malla.X,malla.Y,T, cmap='inferno')
cbar=surf.colorbar(s, shrink=0.5)



plt.show()

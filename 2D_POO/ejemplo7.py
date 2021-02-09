import FVM2D as fvm
import numpy as np
import matplotlib.pyplot as plt

Lx = 10 # meters
Ly = 10
poro = 0.2
TA = 1.00 #
TB = 2.00 # 
TC = 1.00
TD = 5.00
q = 1e+06
rho = 1.0 # kg/m^3
mu = 1.0
k  = 1000 # W/m.K
cT = 1e-4
ux = .1 # m/s
uy = .1 # m/s

Tmax = 0.8
delta_t = 0.0001 # Paso de tiempo
steps = int(Tmax/delta_t)

# Parámetros numéricos
m = 1
Nx = 10 * 2**m + 1
Ny = 10 * 2**m + 1
x = np.linspace(0,100,Nx)
y = np.linspace(0,100,Ny)
dx = Lx / (Nx-1)
dy = Ly / (Ny-1)
t = np.array([0.0005,0.0025,0.0075, 0.0150, 0.0375, 0.8000])
Gamma = k / (poro * mu * cT)

malla = fvm.Mesh2D(Nx,Ny,lengthX=Lx, lengthY=Ly)
print(malla.nodesX,malla.volumesX)
malla.createMesh()
print(malla.dx)

coef = fvm.Coefficients2D(malla.volumesX,malla.volumesY,malla.dx,malla.dy)
coef.alloc()

dif = fvm.Diffusion2D(malla.volumesX,malla.volumesY,malla.dx,malla.dy,Gamma=Gamma)
tem = fvm.Temporal2D(malla.volumesX,malla.volumesY,rho,malla.dx,malla.dy,delta_t)

T = np.ones([coef.nvx,coef.nvy])*0 #Condicion inicial
#Le aplicamos las condiciones de frontera a los coeficientes y a el vector solucion

T[:,0] = TA
T[:,-1] = TD
T[0,:] = TC
T[-1,:] = TB
T_aux=T[1:-1,1:-1].ravel()

# coef.bcDirichlet('LEFT_WALL',TA)
# coef.bcDirichlet('TOP_WALL',TC)
# coef.bcDirichlet('RIGHT_WALL',TD)
# coef.bcDirichlet('BOTTOM_WALL',TB)

for i in range(1, steps+1):
    time_k = i * delta_t
    coef.cleanCoefficients()
    dif.calcCoef()
    tem.calcCoef(T)

    coef.bcDirichlet('LEFT_WALL',TA)
    coef.bcDirichlet('TOP_WALL',TC)
    coef.bcDirichlet('RIGHT_WALL',TD)
    coef.bcDirichlet('BOTTOM_WALL',TB)

    A = fvm.Matrix2D(malla.volumesX,malla.volumesY)
    A.build(coef)

    T_aux = np.linalg.solve(A.A,coef.Su[1:-1,1:-1].flatten())
    T_aux.shape = (Nx-1,Ny-1)
    T[1:-1,1:-1] =T_aux

    
f1 = plt.figure()
c = plt.contourf(malla.X,malla.Y,T,8, alpha=.75,cmap='inferno')
f1.colorbar(c, shrink=1.0)


surf=plt.figure(figsize=(5,4)) 
ax=surf.gca(projection='3d')
s=ax.plot_surface(malla.X,malla.Y,T, cmap='inferno')
cbar=surf.colorbar(s, shrink=0.5)



plt.show()









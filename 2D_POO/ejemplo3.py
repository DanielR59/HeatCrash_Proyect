
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg.linalg import det
nvx = 5
A=np.zeros([nvx,nvx])
U= np.zeros([nvx-1,nvx-1])

for i in range(1,nvx-1):
    for j in range(1,nvx-1):
            U[i-1,j] = 1

print(U[3])


plt.figure()
plt.plot(U)


def hola(ds,dsd,dsds):


    if ds == 1:
        pass







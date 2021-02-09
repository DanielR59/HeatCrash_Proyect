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
k  = 1000 # W/m.K #Gamma
ux = .1 # m/s
uy = .1 # m/s

delta_t = 0.002 # Paso de tiempo
steps = 500


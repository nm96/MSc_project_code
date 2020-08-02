import numpy as np
import scipy.integrate as spi
import solutions as sln
import models as mdl
import pickle

import time; t0 = time.time()

# Integration limits:
T = 2**14
tspan = (0,T)  
N = T*2**4
tt = np.linspace(*tspan,N)
X0 = [0.01,0,0,0]

# Parameter values and model:
eps = 0.2 # Rotor eccentricity
k = 10 # Stiffness coefficient
m = 10 # Mass (per unit length)
h = 1 # Gap width

Om = 4.1
c = 1.0
B = 2.7

model = (mdl.NHSommerfeld2,(Om,h,B))
params = (eps,Om,m,c,k,model)

# Solve:
sol = spi.solve_ivp(sln.dXdt,tspan,X0,t_eval=tt,args=params,method='Radau')

# Write parameter-solution tuple to a binary file:
sol_data = (params, sol)
with open("sol_data.pkl", 'wb') as f:
    pickle.dump(sol_data,f,pickle.HIGHEST_PROTOCOL)

# Save a copy of the file in an external data directory:
fname = "../data/raw/B{:.2f}_Om{:.2f}_c{:.2f}_T{:.0e}_sol_data.pkl".format(B,Om,c,T)
with open(fname, 'wb') as f:
    pickle.dump(sol_data,f,pickle.HIGHEST_PROTOCOL)


tf = time.time()
print("Time taken to compute solution = {:.2f}s".format(tf-t0))

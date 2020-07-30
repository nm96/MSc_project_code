import numpy as np
from scipy.integrate import solve_ivp
from solutions import *
from models import *
import pickle

import time; t0 = time.time()

# Integration limits:
tspan = (0,2**14)  
N = tspan[1]*2**4
tt = np.linspace(*tspan,N)
X0 = [0.01,0,0,0]

# Parameter values and model:
eps = 0.2 # Rotor eccentricity
Om = 4.1 # Driving frequency
m = 10 # Mass (per unit length)
c = 1.5 # Damping coefficient
k = 10 # Stiffness coefficient
h = 1 # Gap width

B = 1.65

model = (NHSommerfeld2,(Om,h,B))
params = (eps,Om,m,c,k,model)

# Solve:
sol = solve_ivp(dXdt,tspan,X0,t_eval=tt,args=params,method='Radau')

# Write parameter-solution tuple to a binary file:
sol_data = (params, sol)
with open("sol_data.pkl", 'wb') as f:
    pickle.dump(sol_data,f,pickle.HIGHEST_PROTOCOL)


tf = time.time()
print("Time taken to compute solution = {:.2f}s".format(tf-t0))

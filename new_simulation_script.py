# A script for integrating the Jeffcott equations using scipy's solve_ivp.

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Rename basic functions and constants for clarity:
cos = np.cos
sin = np.sin
tanh = np.tanh
pi = np.pi

def simple(X):
    return (0,0)

def dXdt(t,X,eps,Om,m,c,k):
    x,dx,y,dy = X
    F_x, F_y = 0,0
    return [dx,
            (F_x + eps*m*cos(Om*t)*Om**2 - c*dx - k*x)/m,
            dy,
            (F_y + eps*m*sin(Om*t)*Om**2 - c*dy - k*y)/m]

eps = 0.03
Om = 4.1
m = 10
c = 0.1
k = 10

params = (eps,Om,m,c,k)
    
tspan = (0,10)    
tt = np.linspace(*tspan,1000)
X0 = [0.1,0.2,0.1,0.2]

sol = solve_ivp(dXdt,tspan,X0,t_eval=tt,args=params)

# Plot the trajectory of (x,y)
plt.plot(sol.y[0],sol.y[2])
plt.show()


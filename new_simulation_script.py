# A script for integrating the Jeffcott equations using scipy's solve_ivp.

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Rename basic functions and constants for clarity
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


# Define parameter values
eps = 0.03
Om = 4.1
m = 10
c = 0.1
k = 10

params = (eps,Om,m,c,k)
    
# Integrat the ODE system over a given time span with given initial conditions
tspan = (0,100)    
tt = np.linspace(*tspan,1000)
X0 = [0.1,0,0,0]

sol = solve_ivp(dXdt,tspan,X0,t_eval=tt,args=params)
solx = sol.y[0]
soly = sol.y[2]

# Plot the trajectory of (x,y) in the stationary frame
fig = plt.figure(1)
ax = fig.add_axes([.1,.1,.8,.8])
ax.plot(solx,soly)

# Plot the trajectory of (x,y) in a rotating frame
rotsolx = cos(Om*tt)*solx - sin(Om*tt)*soly
rotsoly = sin(Om*tt)*solx + cos(Om*tt)*soly

fig = plt.figure(2)
ax = fig.add_axes([.1,.1,.8,.8])
ax.plot(rotsolx,rotsoly)

plt.show()

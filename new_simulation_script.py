# A script for integrating the Jeffcott equations using scipy's solve_ivp.

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

fn = 0  # Initialize figure number for plotting

# Rename basic functions and constants for clarity
cos = np.cos
sin = np.sin
tanh = np.tanh
pi = np.pi

def simple(X,dummy):
    """Model function for the trivial no-forcing model. Dummy parameter
    included in order to fit the convention for model-passing.
    """
    return (0,0)

def dXdt(t,X,eps,Om,m,c,k,h,model):
    # Unpack the components of X:
    x,dx,y,dy = X
    # Define r = distance from stator centre to rotor centre:
    r = (x*x + y*y + eps*eps - 2*eps*(x*cos(Om*t) + y*sin(Om*t)))**0.5
    cPsi = (x - eps*cos(Om*t))/r     # = cos(psi)
    sPsi = (y - eps*sin(Om*t))/r     # = sin(psi)
    # Calculate radial and tangential forces using the given model:
    fn, ft = model[0](X,*model[1])
    # Convert these into forces in the x and y directions:
    F_x = -fn*cPsi + ft*sPsi
    F_y = -fn*sPsi - ft*cPsi
    # Final result (1st order Jeffcott eqns):
    return [dx,
            (F_x + eps*m*cos(Om*t)*Om**2 - c*dx - k*x)/m,
            dy,
            (F_y + eps*m*sin(Om*t)*Om**2 - c*dy - k*y)/m]


# Define parameter values
eps = 0.03  # Rotor eccentricity
Om = 4.1    # Driving frequency
m = 10      # Mass (per unit length)
c = 0.1     # Damping coefficient
k = 10      # Stiffness coefficient
h = 0.03    # Gap width

model = (simple,(0,)) # Model-tuple for the trivial no-force model

params = (eps,Om,m,c,k,h,model) # Package parameters into a tuple
    
# Integrate the ODE system over a given time span with given initial conditions
tspan = (0,100)    
tt = np.linspace(*tspan,1000)
X0 = [0.1,0,0,0]

sol = solve_ivp(dXdt,tspan,X0,t_eval=tt,args=params)
solx = sol.y[0]
soly = sol.y[2]

# Plot the trajectory of (x,y) in the stationary frame
fn += 1; fig = plt.figure(fn)
ax = fig.add_axes([.1,.1,.8,.8])
ax.plot(solx,soly)

# Plot the trajectory of (x,y) in a rotating frame
rotsolx = cos(Om*tt)*solx - sin(Om*tt)*soly
rotsoly = sin(Om*tt)*solx + cos(Om*tt)*soly

fn += 1; fig = plt.figure(fn)
ax = fig.add_axes([.1,.1,.8,.8])
ax.plot(rotsolx,rotsoly)

# Integration with VdH model

def VdH(X,h,k_c):
    """Model function for the Van der Heijden model. Takes two parameters h and
    k_c, with the latter being model-specific"""
    x = X[0]
    y = X[2]
    r_H = (x*x + y*y)**0.5
    delta_H = (r_H-h)*(tanh(1000*(r_H-h))+1)/2
    return (k_c*delta_H,0)

# Define the model
k_c = 100 # Required model-specific parameter
model = (VdH,(h,k_c)) # This is the standard form for a model;
# A tuple containing a model function with the form f(X,mparams) and a tuple
# mparams of model parameters.

params = (eps,Om,m,c,k,h,model)
    
# Integrate the ODE system over a given time span with given initial conditions
tspan = (0,100)    
tt = np.linspace(*tspan,1000)
X0 = [0.1,0,0,0]

sol = solve_ivp(dXdt,tspan,X0,t_eval=tt,args=params)
solx = sol.y[0]
soly = sol.y[2]

# Plot the trajectory of (x,y) in the stationary frame
fn += 1; fig = plt.figure(fn)
ax = fig.add_axes([.1,.1,.8,.8])
ax.plot(solx,soly)

# Plot the trajectory of (x,y) in a rotating frame
rotsolx = cos(Om*tt)*solx - sin(Om*tt)*soly
rotsoly = sin(Om*tt)*solx + cos(Om*tt)*soly

fn += 1; fig = plt.figure(fn)
ax = fig.add_axes([.1,.1,.8,.8])
ax.plot(rotsolx,rotsoly)

plt.show()

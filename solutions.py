"""Module containing functions for solving the Jeffcott rotordynamic equations,
processing and plotting the resulting solutions. Specific nonlinear models to
be incorporated in the Jeffcott equations are in the the 'models' module.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.fft import rfft
import matplotlib.pyplot as plt

# Rename basic functions and constants for clarity
cos = np.cos
sin = np.sin
tanh = np.tanh
pi = np.pi

def dXdt(t,X,eps,Omega,m,c,k,h,model):
    """Right hand side of the Jeffcott equations in first order form, to be
    passed to the numerical integrator 'solve_ivp'.
    """
    # Unpack the components of X:
    x,dx,y,dy = X
    # Define r = distance from stator centre to rotor centre:
    r = (x*x + y*y + eps*eps - 2*eps*(x*np.cos(Omega*t) + y*np.sin(Omega*t)))**0.5
    cPsi = (x - eps*np.cos(Omega*t))/r     # = np.cos(psi)
    sPsi = (y - eps*np.sin(Omega*t))/r     # = np.sin(psi)
    # Calculate radial and tangential forces using the given model:
    fn, ft = model[0](X,*model[1])
    # Convert these into forces in the x and y directions:
    F_x = -fn*cPsi + ft*sPsi
    F_y = -fn*sPsi - ft*cPsi
    # Final result (1st order Jeffcott eqns):
    return [dx, (F_x + eps*m*np.cos(Omega*t)*Omega**2 - c*dx - k*x)/m, dy, (F_y
        + eps*m*np.sin(Omega*t)*Omega**2 - c*dy - k*y)/m]


def solxy(sol):
    """Return an array of solution x and y values from a solution object 'sol'
    produced by solve_ivp
    """
    return sol.y[(0,2),:]

def rotsolxy(sol,Omega):
    """Return an array of x and y solution values in a frame rotating at angular
    velocity Omega. Requires a solution object 'sol' output by solve_ivp as
    well as the relevant angular velocity.
    """
    solx = sol.y[0]
    soly = sol.y[2]
    tt = sol.t
    return np.array([np.cos(Omega*tt)*solx - np.sin(Omega*tt)*soly,
    np.sin(Omega*tt)*solx + np.cos(Omega*tt)*soly])

def mod2dfft(sol):
    solx = sol.y[0]
    soly = sol.y[2]
    fftx = rfft(solx)
    ffty = rfft(soly)
    return (np.abs(fftx)**2 + np.abs(ffty)**2)**0.5


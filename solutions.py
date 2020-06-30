"""Module containing functions for solving the Jeffcott rotordynamic equations,
processing and plotting the resulting solutions. Specific nonlinear models to
be incorporated in the Jeffcott equations are in the the 'models' module.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.fft import fft
import matplotlib.pyplot as plt

# Rename basic functions and constants for clarity
cos = np.cos
sin = np.sin
tanh = np.tanh
pi = np.pi

def dXdt(t,X,eps,Om,m,c,k,h,model):
    """Right hand side of the Jeffcott equations in first order form, to be
    passed to the numerical integrator 'solve_ivp'.
    """
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


def sol_xy(sol):
    return sol.y[(0,2),:]

def rotsol_xy(sol,Om):
    solx = sol.y[0]
    soly = sol.y[2]
    tt = sol.t
    return np.array([cos(Om*tt)*solx - sin(Om*tt)*soly,
    sin(Om*tt)*solx + cos(Om*tt)*soly])

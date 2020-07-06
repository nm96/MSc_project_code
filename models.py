"""Module containing functions which encode various nonlinear forcing models
which can be incorporated into the Jeffcott equations.

Each model function takes the phase-space position vector X along with some
combination of parameters and returns the tuple 
(fn,ft) = (normal force, tangential force).
"""

import numpy as np

# Rename basic functions and constants for clarity
cos = np.cos
sin = np.sin
tanh = np.tanh
pi = np.pi

def simple(X,dmy):
    """Model function for the trivial no-forcing model. Dummy parameter
    included in order to fit the convention for model-passing.
    Returns the tuple (fn,ft) = (normal force, tangential force).
    """
    return (0,0)

def VdH(X,h,k_c):
    """Model function for the Van der Heijden model. Takes two parameters h and
    k_c, with the latter being model-specific
    Returns the tuple (fn,ft) = (normal force, tangential force).
    """
    x = X[0]
    y = X[2]
    r = (x*x + y*y)**0.5
    #delta = (r-h)*(tanh(1000*(r-h))+1)/2
    delta = (r-h)*np.heaviside(r-h,0.5)
    return (k_c*delta,0)

def NHSommerfeld(X,Omega,h,mu,b,R2):
    """Model function for the Naive Half-Sommerfeld lubrication model. See
    preliminary report 3.4.3 for details. Takes global parameters Omega and h
    along with model-specific parameters mu,b,R2.
    """
    pass

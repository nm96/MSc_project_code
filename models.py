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
    """
    return (0,0)

def VdH(X,h,k_c):
    """Model function for the Van der Heijden model. Takes two parameters h and
    k_c, with the latter being model-specific
    """
    x = X[0]
    y = X[2]
    r = (x*x + y*y)**0.5
    delta = (r-h)*np.heaviside(r-h,0.5)
    return (k_c*delta,0)

def Hua(X,h,k_c,alpha,mu):
    """Model function for the Hua model. Takes three parameters: h, k_c, alpha,
    mu.
    REQUIRES A DIRECTION-DEPENDENT FRICTION FORCE!
    """
    x = X[0]
    y = X[2]
    r = (x*x + y*y)**0.5
    delta = (r-h)*np.heaviside(r-h,0.5)
    return (k_c*delta + alpha*delta**2,mu*(k_c*delta + alpha*delta**2))

def NHSommerfeld(X,Omega,h,mu,b,R2):
    """Model function for the Naive Half-Sommerfeld lubrication model. See
    preliminary report 3.4.3 for details. Takes global parameters Omega and h
    along with model-specific parameters mu,b,R2.
    """
    x,dx,y,dy = X
    r = (x*x + y*y)**0.5
    dr = (x*dx + y*dy)/r
    a = r/h
    da = dr/h
    dpsi = (x*dy - y*dx)/r**2
    fn = -(12*mu*b*R2**3/h**2)*(a**2*(Omega-2*dpsi)/((2+a**2)*(1-a**2)) +
            da*(1-a**2)**(-1.5)*(pi/2 - 8/(2+a**2))
    ft = (12*mu*b*R2**3/h**2)*(pi*a*(Omega-2*dpsi)/(2*(2+a**2)(1-a**2)**0.5) +
        2*a*da/((2+a**2)*(1-a**2)))

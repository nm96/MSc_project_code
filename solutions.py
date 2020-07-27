"""Module containing functions for solving the Jeffcott rotordynamic equations,
processing and plotting the resulting solutions. Specific nonlinear models to
be incorporated in the Jeffcott equations are in the the 'models' module.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.fft import fft
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
    NEW VERSION: The Jeffcott equations should describe the motion of the
    rotor's geometric center rather than its center of mass.
    """
    # Unpack the components of X:
    x,dx,y,dy = X
    # Define r = distance from stator centre to rotor centre:
    r = (x*x + y*y)**0.5
    cPsi = x/r    # = np.cos(psi)
    sPsi = y/r    # = np.sin(psi)
    # Calculate radial and tangential forces using the given model:
    fn, ft = model[0](X,*model[1])
    # Convert these into forces in the x and y directions:
    F_x = -fn*cPsi + ft*sPsi
    F_y = -fn*sPsi - ft*cPsi
    # Final result (1st order Jeffcott eqns):
    return [dx, (F_x + eps*m*np.cos(Omega*t)*Omega**2 - c*dx - k*x)/m, dy, (F_y
        + eps*m*np.sin(Omega*t)*Omega**2 - c*dy - k*y)/m]

def dXdt_old(t,X,eps,Omega,m,c,k,h,model):
    """Right hand side of the Jeffcott equations in first order form, to be
    passed to the numerical integrator 'solve_ivp'.
    OLD (INCORRECT) VERSION
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
    produced by solve_ivp. Automatically removes transients (first half of the
    time range)
    """
    N = len(sol.t)
    return sol.y[(0,2),N//2:]

def rotsolxy(sol,Omega):
    """Return an array of x and y solution values in a frame rotating at
    angular velocity Omega. Requires a solution object 'sol' output by
    solve_ivp as well as the relevant angular velocity. Automatically removes
    transients (first half of the time range).
    """
    N = len(sol.t)
    solx = sol.y[0,N//2:]
    soly = sol.y[2,N//2:]
    tt = sol.t[N//2:]
    return np.array([np.cos(Omega*tt)*solx + np.sin(Omega*tt)*soly,
        - np.sin(Omega*tt)*solx + np.cos(Omega*tt)*soly])

def solxy2(sol):
    """Return an array of solution x and y values from a solution object 'sol'
    produced by solve_ivp. This version keeps transients.
    """
    N = len(sol.t)
    return sol.y[(0,2),:]

def rotsolxy2(sol,Omega):
    """Return an array of x and y solution values in a frame rotating at
    angular velocity Omega. Requires a solution object 'sol' output by
    solve_ivp as well as the relevant angular velocity. This version keeps
    transients.
    """
    N = len(sol.t)
    solx = sol.y[0]
    soly = sol.y[2]
    tt = sol.t
    return np.array([np.cos(Omega*tt)*solx + np.sin(Omega*tt)*soly,
        - np.sin(Omega*tt)*solx + np.cos(Omega*tt)*soly])

def transformed(sol):
    """Calculate the logarithm of the combined absolute values of the fourier
    transforms of a solution's x and y components. Returns both the transform
    result and the corresponding vector of frequencies.  
    """
    N = len(sol.t)
    max_freq = min(10,1/sol.t[1])
    freq_scale = np.arange(0,max_freq,4*pi/sol.t[-1])
    solx = sol.y[0]
    soly = sol.y[2]
    fftx = rfft(solx[N//2:])[:len(freq_scale)]
    ffty = rfft(soly[N//2:])[:len(freq_scale)]
    return (freq_scale, np.log((np.abs(fftx)**2 + np.abs(ffty)**2)**0.5))

def transformed2(sol):
    """Calculate the logarithm of the combined absolute values of the fourier
    transforms of a solution's x and y components. Returns both the transform
    result and the corresponding vector of frequencies.  
    This version uses a hanning window to smooth the solution data.
    """
    N = len(sol.t)
    max_freq = min(10,1/sol.t[1])
    freq_scale = np.arange(0,max_freq,4*pi/sol.t[-1])
    solx = sol.y[0][N//2:]
    soly = sol.y[2][N//2:]
    fftx = rfft(solx*np.hanning(len(solx)))[:len(freq_scale)]
    ffty = rfft(soly*np.hanning(len(soly)))[:len(freq_scale)]
    return (freq_scale, np.log((np.abs(fftx)**2 + np.abs(ffty)**2)**0.5))

def PSD(sol):
    """New power spectrum density function. Essentially equivalent to the old
    transform functions above but using only one component (x(t)) of the
    solution and applying a Hanning window before taking the Fourier
    transform.
    """
    N0 = len(sol.t)
    x = sol.y[0][N0//2:]
    N = len(x)
    T = sol.t[-1]
    w = np.hanning(N)
    P = (2/N)*abs(fft(w*x))
    om = np.arange(N)*4*pi/T
    return (om,P)

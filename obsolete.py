# Old rarely used functions - stashing them here rather than using git
# properly!

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


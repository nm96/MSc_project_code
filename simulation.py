import numpy as np
import scipy.integrate as spi
import matplotlib
import matplotlib.pyplot as plt
import pickle
import time
import sys

# Set a nice font for the plots when on Linux
if sys.platform == 'linux':
    matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    matplotlib.rc('text', usetex=True)



class Simulation:


    # Default physical parameters:
    eps = 0.2 
    k = 10
    m = 10
    h = 1
    Om = 4.1
    c = 0.5
    B = 0.5

    # Default integration parameters
    X0 = [0.01,0,0,0]
    T = 2**14
    N = T*2**4
    rtol = 1e-4
    atol = 1e-8
    
    def model(self,X):
        """Model function for the Naive Half-Sommerfeld lubrication model. See
        preliminary report 3.4.3 for details.
        """
        x,dx,y,dy = X
        r = (x*x + y*y)**0.5
        dr = (x*dx + y*dy)/r
        a = r/self.h
        da = dr/self.h
        dpsi = (x*dy - y*dx)/r**2
        fn = (-(self.B/self.h**2)*(a**2*(self.Om-2*dpsi)/
            ((2+a**2)*(1-a**2))+ da*(1-a**2)**(-1.5)*(np.pi/2 - 8/(2+a**2)))) 
        ft = ((self.B/self.h**2)*(np.pi*a*(self.Om-2*dpsi)/
            (2*(2+a**2)*(1-a**2)**0.5) + 2*a*da/((2+a**2)*(1-a**2))))
        return (fn,ft)

    def dXdt(self,t,X):
        """Right hand side of the Jeffcott equations in first order form, to be
        passed to the numerical integrator 'solve_ivp'. 
        """
        # Unpack the components of X:
        x,dx,y,dy = X
        # Define r = distance from stator centre to rotor centre:
        r = (x*x + y*y)**0.5
        cPsi = x/r    # = np.cos(psi)
        sPsi = y/r    # = np.sin(psi)
        # Calculate radial and tangential forces using the given model:
        fn, ft = self.model(X)
        # Convert these into forces in the x and y directions:
        F_x = -fn*cPsi + ft*sPsi
        F_y = -fn*sPsi - ft*cPsi
        # Final result (1st order Jeffcott eqns):
        return [dx, (F_x + self.eps*self.m*np.cos(self.Om*t)*self.Om**2 -
            self.c*dx - self.k*x)/self.m, dy, (F_y +
                self.eps*self.m*np.sin(self.Om*t)*self.Om**2 - self.c*dy -
                self.k*y)/self.m]
        
    def solve(self):
        """Method for actually integrating the equations."""
        t0 = time.time()
        X0 = self.X0
        T = self.T
        N = self.N
        tspan = (0,T)
        sol = spi.solve_ivp(self.dXdt,tspan,X0,t_eval=np.linspace(0,T,N),
                rtol=self.rtol,atol=self.atol,method='Radau')
        self.t = sol.t
        self.X = sol.y
        tf = time.time()
        print("Solution time = {:.2f}s".format(tf-t0))

    def transform(self):
        """Method for producing a power spectrum density"""
        t0 = time.time()
        N1 = self.N//2
        T = self.T
        x = self.X[0][N1:]
        w = np.hanning(N1)
        self.P = (2/N1)*abs(np.fft.fft(w*x))**2
        self.om = np.arange(N1)*4*np.pi/T
        tf = time.time()
        print("Transform time = {:.2f}s".format(tf-t0))

    def psd_plot(self,fn=1,om_max=10):
        fig, ax = plt.subplots(num=fn)
        ax.set_xlim([0,om_max])
        om_nat = (self.k/self.m)**0.5
        ax.axvline(om_nat,ls='--',c='g',label=r"$\omega_{nat}$")
        ax.axvline(self.Om,ls='--',c='r',label=r"$\Omega$")
        ax.semilogy(self.om,self.P,c='k')
        locmaj = matplotlib.ticker.LogLocator(base=100,numticks=30) 
        ax.yaxis.set_major_locator(locmaj)
        locmin = matplotlib.ticker.LogLocator(base=100,subs=(0.2,0.4,0.6,0.8),
                numticks=50)
        ax.yaxis.set_minor_locator(locmin)
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.grid()
        ts = r"$\beta={:.2f}$, $\Omega={:.2f}$, $c={:.2f}$"
        ts = ts.format(self.B,self.Om,self.c)
        ax.set_title(ts)
        ax.set_ylabel("$P(\omega)$",rotation=0)
        ax.yaxis.labelpad = 20
        ax.set_xlabel("$\omega$")
        ax.legend()
        plt.tight_layout()

    def phase_plot(self,fn=1):
        N1 = self.N//2
        x = self.X[0,N1:]
        y = self.X[2,N1:]
        tt = self.t[N1:]
        arr = np.array([np.cos(self.Om*tt)*x + np.sin(self.Om*tt)*y,
        - np.sin(self.Om*tt)*x + np.cos(self.Om*tt)*y])
        fig, ax = plt.subplots(num=fn)
        ax.plot(*arr)
        ax.set_aspect('equal')
        ts = r"$\beta={:.2f}$, $\Omega={:.2f}$, $c={:.2f}$"
        ts = ts.format(self.B,self.Om,self.c)
        ax.set_title(ts)

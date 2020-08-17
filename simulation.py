import numpy as np
import scipy.integrate as spi
import scipy.signal as sps
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
    #N = T*2**4
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
        N = T*2**6
        self.N = N
        tspan = (0,T)
        sol = spi.solve_ivp(self.dXdt,tspan,X0,t_eval=np.linspace(0,T,N),
                rtol=self.rtol,atol=self.atol,method='Radau')
        self.t = sol.t
        self.X = sol.y
        tf = time.time()
        print("Solution time = {:.2f}s".format(tf-t0))

    def transform(self,R=False):
        """Method for producing a power spectrum density"""
        t0 = time.time()
        N = len(self.t)
        N1 = N//2
        T = self.T
        if R:
            self.rotate()
            x = self.RX[0][N1:]
            w = np.hanning(N1)
            self.P = (2/N1)*abs(np.fft.fft(w*x))**2
            self.om = np.arange(N1)*4*np.pi/T
        else:
            x = self.X[0][N1:]
            w = np.hanning(N1)
            self.P = (2/N1)*abs(np.fft.fft(w*x))**2
            self.om = np.arange(N1)*4*np.pi/T
        tf = time.time()
        print("Transform time = {:.2f}s".format(tf-t0))

    def psd_plot(self,fn=1,om_max=10,R=False,fs=[6,4]):
        self.transform(R)
        fig, ax = plt.subplots(num=fn,figsize=fs)
        ax.set_xlim([0,om_max])
        om_nat = (self.k/self.m)**0.5
        ax.semilogy(self.om,self.P,c='k',linewidth=1)
        locmaj = matplotlib.ticker.LogLocator(base=100,numticks=30) 
        ax.yaxis.set_major_locator(locmaj)
        locmin = matplotlib.ticker.LogLocator(base=100,subs=(0.2,0.4,0.6,0.8),
                numticks=50)
        ax.yaxis.set_minor_locator(locmin)
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.grid()
        ax.set_ylabel("$P(\omega)$",rotation=0)
        ax.yaxis.labelpad = 20
        ax.set_xlabel("$\omega$")
        ts = r"$\beta={:.2f}$, $\Omega={:.2f}$, $c={:.2f}$"
        ts = ts.format(self.B,self.Om,self.c)
        ax.set_title(ts)
        if R == False:
            ax.axvline(om_nat,ls='--',c='g',label=r"$\omega_{nat}$")
            ax.axvline(self.Om,ls='--',c='r',label=r"$\Omega$")
            ax.legend()
        plt.tight_layout()
        fnm = "../plots/B{:.2f}Om{:.2f}c{:.2f}psdplot.eps"
        self.psd_plot_filename = fnm.format(self.B,self.Om,self.c)

    def find_peaks(self,om_max=10):
        P = self.P[self.om < om_max]
        om = self.om[self.om < om_max]
        self.peaks = om[sps.find_peaks(P,prominence=1e-5)[0]]

    def rotate(self):
        """Transform all components of the solution to the rotating frame"""
        x, dx, y, dy = self.X
        t = self.t
        Om = self.Om
        Rx = x*np.cos(Om*t) + y*np.sin(Om*t)
        Ry = -x*np.sin(Om*t) + y*np.cos(Om*t)
        Rdx = dx*np.cos(Om*t) + dy*np.sin(Om*t) + Om*Ry
        Rdy = -dx*np.sin(Om*t) + dy*np.cos(Om*t) - Om*Rx
        self.RX = np.array([Rx,Rdx,Ry,Rdy])

    def phase_plot(self,fn=1,d=2,fs=[6,4]):
        self.rotate()
        self.N = len(self.t)
        N1 = self.N - self.N//d
        Rx,Rdx,Ry,Rdy = self.RX[:,N1:]
        fig = plt.figure(fn,figsize=fs)
        ax = fig.add_subplot(121)
        ax.plot(Rx,Ry,'.',ms=0.05)
        ax.set_aspect('equal')
        ax.set_xlabel(r"\Large $\tilde{x}$")
        ax.set_ylabel(r"\Large $\tilde{y}$",rotation=0)
        ax = fig.add_subplot(122)
        ax.plot(Rdx,Rdy,'.',ms=0.05)
        ax.set_aspect('equal')
        ax.set_xlabel(r"\Large $\dot{\tilde{x}}$")
        ax.set_ylabel(r"\Large $\dot{\tilde{y}}$",rotation=0)
        ts = r"""Solution trajectories in the rotating frame 
        $\beta={:.2f}$, $\Omega={:.2f}$, $c={:.2f}$"""
        ts = ts.format(self.B,self.Om,self.c)
        fig.suptitle(ts)
        plt.tight_layout()
        fnm = "../plots/B{:.2f}Om{:.2f}c{:.2f}phaseplot.eps"
        self.phase_plot_filename = fnm.format(self.B,self.Om,self.c)

    def first_return_plot(self):
        t0 = time.time()
        X0 = self.X0
        T = self.T
        tspan = (0,T)
        T_Om = 2*np.pi/self.Om
        t_eval = np.arange(T//T_Om)*T_Om
        sol = spi.solve_ivp(self.dXdt,tspan,X0,t_eval=t_eval,
                rtol=self.rtol,atol=self.atol,method='Radau')
        tf = time.time()
        print("Solution time = {:.2f}s".format(tf-t0))
        fig, ax = plt.subplots(num=1)
        ax.scatter(sol.y[0],sol.y[2],marker='x',c=list(t_eval))
        ax.set_aspect('equal')
        ts = r"$\beta={:.2f}$, $\Omega={:.2f}$, $c={:.2f}$"
        ts = ts.format(self.B,self.Om,self.c)
        ax.set_title(ts)
        fig, ax = plt.subplots(num=2)
        ax.scatter(sol.y[1],sol.y[3],marker='x',c=list(t_eval))
        ax.set_aspect('equal')
        ts = r"$\beta={:.2f}$, $\Omega={:.2f}$, $c={:.2f}$"
        ts = ts.format(self.B,self.Om,self.c)
        ax.set_title(ts)


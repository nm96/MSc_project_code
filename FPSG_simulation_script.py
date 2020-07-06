# Study group code converted into a script - to be used as a basis for later
# code, it does not actually run properly.
# --------------------------------------------------------------------------

# coding: utf-8

# # Mathematical models
# 
# Jupyter notebook code implementing different models of rotor dynamics.

#  

# In[15]:


# Import of necessary libraries and simplifying notation
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from scipy.integrate import ode
import matplotlib.pyplot as plt
import math
import numpy as np
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib import colors as mcolors

cos = math.cos
sin = math.sin
tanh = math.tanh
pi = math.pi


# ## Defining system for solvers

# In[16]:


# Defining system for ode_int (solver with fixed size time step)
def system(X,t):
    global model
    model = model.upper() # just in case 
    x,dx,y,dy = X
    
    r = (x*x + y*y + eps*eps - 2*eps*(x*cos(Omega*t) + y*sin(Omega*t)))**0.5
    
    # 
    delta = (r - h)*(tanh(1000*(r-h))+1)/2  # or use r - h if r>h else 0  for sharper transition
    
    # calculating derivative of r (r_dot)
    B = x*cos(Omega*t) + y*sin(Omega*t)
    A = x*x + y*y + eps*eps - 2*eps*B
    B_dot = dx*cos(Omega*t) + dy*sin(Omega*t) - Omega*(x*sin(Omega*t) + y*cos(Omega*t))
    A_dot = 2*x*dx + 2*y*dy - 2*eps*B_dot  
    r_dot = 0.5*(A_dot)*(A**-0.5) 
    
    delta_dot = r_dot if r>h else 0
    
      
    cPsi = (x - eps*cos(Omega*t))/r
    sPsi = (y - eps*sin(Omega*t))/r
    
    if model == 'S':
        # simple - no additional force
        fn = 0
        ft = 0
        
    elif model == 'G':
        # Groll
        fn = k_c * (delta**1.5) * (1 + 1.5 * c_c * delta_dot)
        ft = mu * fn
    
    elif model == 'SG':
        # Simple Groll
        fn = k_c * delta + c_c * delta_dot
        ft = mu * fn
    
    elif model == 'VDH':
        # Van der Heijden - slightly different r and delta (why??)
        r_H = (x*x + y*y)**0.5
        delta_H = (r_H - h)*(tanh(1000*(r_H-h))+1)/2
        fn = k_c * delta_H
        ft = 0
        
    #elif model == 'MM':
    # modified Muszynska - SKIP NOW
    #     fn =
    #     ft = mu * fn
    
    elif model == 'OM':
        # original Muszynkska
        fn = k_c * r if ((Omega*t)%(2*pi) >= omega) else 0
        ft = mu * fn
    
    elif model == 'H':
        # Hua
        fn = k_c * delta + alpha * delta**2
        ft = mu * fn
    
    elif model == 'L':
        # Lubrication
        Psi_dot = ((dy + eps*(Omega**2)*sin(Omega*t))*r - (y - eps*Omega*cos(Omega*t))*r_dot)/((r**2) * cPsi)
        fn = - (2*pi*eta*(R2**3)*r_dot)/((h**2-r**2)**1.5)
        ft = (12*pi*eta*(R2**3)*(Omega - 2*Psi_dot)*r)/((h**2 - r**2)**0.5 * (2*h**2 + r**2))
    
    else:
        print("Wrong choice of model. Continuing without any additional forces.")
        fn = 0
        ft = 0
    
    F_x = -fn*cPsi + ft*sPsi
    F_y = -fn*sPsi - ft*cPsi
    
    dXdt = [dx,
            (F_x + eps*m*cos(Omega*t)*(Omega**2) + c_1*m*Omega*dy - c_2*dx - k_1*x - k_2*y)/m,
            dy,
            (F_y + eps*m*sin(Omega*t)*(Omega**2) - c_1*m*Omega*dx - c_2*dy + k_2*x - k_1*y)/m           
            ]
    return dXdt

def system_ivp(t,X):
    global model
    model = model.upper() # just in case 
    x,dx,y,dy = X
    
    r = (x*x + y*y + eps*eps - 2*eps*(x*cos(Omega*t) + y*sin(Omega*t)))**0.5
    
    # 
    delta = (r - h)*(tanh(1000*(r-h))+1)/2  # or use r - h if r>h else 0  for sharper transition
    
    # calculating derivative of r (r_dot)
    B = x*cos(Omega*t) + y*sin(Omega*t)
    A = x*x + y*y + eps*eps - 2*eps*B
    B_dot = dx*cos(Omega*t) + dy*sin(Omega*t) - Omega*(x*sin(Omega*t) + y*cos(Omega*t))
    A_dot = 2*x*dx + 2*y*dy - 2*eps*B_dot  
    r_dot = 0.5*(A_dot)*(A**-0.5) 
    
    delta_dot = r_dot if r>h else 0
    
    cPsi = (x - eps*cos(Omega*t))/r     # = cos(psi)
    sPsi = (y - eps*sin(Omega*t))/r     # = sin(psi)
    
    if model == 'S':
        # simple - no additional force
        fn = 0
        ft = 0
        
    elif model == 'G':
        # Groll
        fn = k_c * (delta**1.5) * (1 + 1.5 * c_c * delta_dot)
        ft = mu * fn
    
    elif model == 'SG':
        # Simple Groll
        fn = k_c * delta + c_c * delta_dot
        ft = mu * fn
    
    elif model == 'VDH':
        # Van der Heijden - slightly different r and delta
        r_H = (x*x + y*y)**0.5
        delta_H = (r_H - h)*(tanh(1000*(r_H-h))+1)/2
        fn = k_c * delta_H
        ft = 0
        
    #elif model == 'MM':
    # modified Muszynska - SKIP NOW
    #     fn =
    #     ft = mu * fn
    
    elif model == 'OM':
        # original Muszynkska
        fn = k_c * r if ((Omega*t)%(2*pi) >= omega) else 0
        ft = mu * fn
    
    elif model == 'H':
        # Hua
        fn = k_c * delta + alpha * delta**2
        ft = mu * fn
    
    elif model == 'L':
        # Lubrication
        Psi_dot = ((dy + eps*(Omega**2)*sin(Omega*t))*r - (y - eps*Omega*cos(Omega*t))*r_dot)/((r**2) * cPsi)
        epsilon =  10**(-2)
        z = h - r - epsilon
        f = (z*(0.5 + 0.5*tanh(z/epsilon)) + epsilon)*(r+h)
        fn = - (2*pi*eta*(R2**3)*r_dot)/(f**1.5) 
        ft = (12*pi*eta*(R2**3)*(Omega - 2*Psi_dot)*r)/((2*h**2 + r**2)* f**0.5)
    
    else:
        print("Wrong choice of model. Continuing without any additional forces.")
        fn = 0
        ft = 0
    
    F_x = -fn*cPsi + ft*sPsi
    F_y = -fn*sPsi - ft*cPsi
    
    
    dXdt = [dx,
            (F_x + eps*m*cos(Omega*t)*(Omega**2) + c_1*m*Omega*dy - c_2*dx - k_1*x - k_2*y)/m,
            dy,
            (F_y + eps*m*sin(Omega*t)*(Omega**2) - c_1*m*Omega*dx - c_2*dy + k_2*x - k_1*y)/m           
            ]
    return dXdt


# ## Time vector and initial condition

# In[17]:


# time vector
t0 = 1000  # end time
dt = 0.01  # desired spacing
t = np.linspace(0, t0, int(t0/dt+1))
t_len = len(t)

# initial condition
X0 = [0.1,0,0,0]


# ## Choosing model, defining parameters

# In[ ]:


# S for simple, G for Groll, SG for Simple Groll, 
# VDH for van der Heijden, OM for original Muszynska, 
# H for Hua, L for Lubrication
model = 'VDH' 

# Parameters of the system
m = 10 # mass
Omega = 4.1 # driving frequency
eps = 0.03 # distance from center of rotation to center of mass
c_1 = 0 # constant for gyroscopic terms
c_2 = 0.1 # damping
k_1 = 10 # stiffness
k_2 = 0 # stiffness cross coupling effects
h = 0.05 # gap parameter
omega_nat = (k_1/m)**0.5 # natural frequency of the rotor

# Parameters of the different models
mu = 0.1
k_c = 100

c_c = 0.2 #0.1 # for Groll
omega = 10 # original Muzynska
alpha = 500 # Hua
eta = 10**(-5) # Lubrication
R2 = h*0.1 # Lubrication


# ## Actual solving

# In[ ]:


# choosing the solver - comment the appropriate line
#solver = 'odeint'
solver = 'solve_ivp'

#choosing the method for solve_ivp
    # stiff problems
# method = 'Radau'
# method = 'BDF'
    # non-stiff problems
method = 'RK45'
# method = 'RK23'
    # universal Fortran method
#method = 'LSODA'

if solver == 'odeint':
    sol = odeint(system, X0, t)
    
elif solver == 'solve_ivp':
    ivp_sol = solve_ivp(system_ivp,[0,t0], X0, method=method, rtol=1e-6, t_eval=t, dense_output=False, events=None, vectorized=True)
    sol1 = ivp_sol.y
    t1 = ivp_sol.t
    sol = np.zeros([t_len,4])
    for I in [0,1,2,3]:
         sol[:,I] = np.interp(t, t1, sol1[I,:])

else:
    print('Error. I do not recognize required solver.')


# In[ ]:


# defining solutions for easier manipulation
disc_param = 1.2
disc = int(t_len//disc_param) # 'percentage' of data discard as transient period
t_new = t[disc:]
solx = sol[:,0]
solx = solx[disc:] # discarding the transient period
soldx = sol[:,1]
soly = sol[:,2]
soly = soly[disc:] # discarding the transient period
soldy = sol[:,3]      

# rotating frame
rotsol1x = np.cos(Omega*t_new)*solx - np.sin(Omega*t_new)*soly
rotsol1y = np.sin(Omega*t_new)*solx + np.cos(Omega*t_new)*soly
rotsol2x = np.cos(Omega*t_new)*solx + np.sin(Omega*t_new)*soly
rotsol2y = -np.sin(Omega*t_new)*solx + np.cos(Omega*t_new)*soly

# plotting
plt.figure(figsize=(10,20))
plt.subplot(211)
plt.title('Phase diagram')
plt.plot(solx,soly)
plt.xlabel('$x$')
plt.ylabel('$y$')

plt.subplot(212)
plt.title('Phase diagram - rotating frame')
plt.plot(rotsol1x,rotsol1y)
plt.xlabel('$x$')
plt.ylabel('$y$')

#plt.savefig(model+'_phase.png')

plt.show()


# ## Calculating FFT (takes the longest...)

# In[ ]:


fftx = np.fft.fft(solx)
fftx2 = (fftx * np.conjugate(fftx))**0.5
ffty = np.fft.fft(soly)
ffty2 =(ffty * np.conjugate(ffty))**0.5

# auxiliary for correct scale for plotting
end_point = 5
x_scale = np.arange(0,end_point,2*pi/(t_new[-1]-t_new[0]))

# plotting
plt.figure(figsize=(10,20))
plt.subplot(211)
plt.title('FFT semilog plot')
plt.plot(x_scale,fftx2[:len(x_scale)], label = 'FFT x')
plt.plot(x_scale,ffty2[:len(x_scale)],label = 'FFT y')
plt.grid()
plt.plot([Omega,Omega],[0,10**4],'r-')
plt.plot([omega_nat,omega_nat],[0,10**4],'r-')
plt.yscale("log")
plt.legend(loc='best')
plt.text(Omega+0.1,10**1,'$\Omega$',fontsize=20)
plt.text(omega_nat+0.1,10**1,'$\Omega_{nat}$',fontsize=20)

plt.subplot(212)
plt.title('Real and imaginary parts of FFT of $x$ and $y$ coordinates')
plt.plot([Omega,Omega],[0,10**4],'r-')
plt.plot([omega_nat,omega_nat],[0,10**4],'r-')
plt.plot(x_scale,np.real(fftx[:len(x_scale)]),label='Re FFT(x)')
plt.plot(x_scale,np.real(ffty[:len(x_scale)]),label='Re FFT(y)')
plt.plot(x_scale,np.imag(fftx[:len(x_scale)]),label='Im FFT(x)')
plt.plot(x_scale,np.imag(ffty[:len(x_scale)]),label='Im FFT(y)')
plt.grid()
plt.yscale("log")
plt.legend(loc='best')
plt.text(Omega+0.1,10**3,'$\Omega$',fontsize=20)
plt.text(omega_nat+0.1,10**3,'$\Omega_{nat}$',fontsize=20)

#plt.savefig(model+'_fft.png')
plt.show()


# ## Making a waterfall plot
# 
# For a waterfall plot, one needs to consider a values generated in a for loop as indicated bellow (example is fow varying frequencies for van der Heijden model

# In[18]:


# time vector
t0 = 1000  # end time
dt = 0.001  # desired spacing
t = np.linspace(0, t0, int(t0/dt+1))
t_len = len(t)

# initial condition
X0 = [0.1,0.2,0.1,0.2]

# S for simple, G for Groll, SG for Simple Groll, 
# VDH for van der Heijden, OM for original Muszynska, 
# H for Hua, L for Lubrication
model = 'VDH' 

# Parameters of the system
m = 10 # mass
Omega = 1 # driving frequency(?)
eps = 0.03 # distance from center of rotation to center of mass
c_1 = 0 # constant for gyroscopic terms
c_2 = 0.01 # damping
k_1 = 10 # stiffness
k_2 = 0 # stiffness cross coupling effects
h = 0.3 # gap parameter
omega_nat = (k_1/m)**0.5 # natural frequency of the rotor

# Parameters of the different models
mu = 0.1
k_c = 100

c_c = 0.2 #0.1 # for Groll
omega = 10 # original Muzynska
alpha = 500 # Hua
eta = 10**(-5) # Lubrication
R2 = h*0.1 # Lubrication


method = 'RK45'

maxim = 25

# definition of the truncated time is moved here for simplicity
disc_param = 1.2
disc = int(t_len//disc_param) # 'percentage' of data discard as transient period
t_new = t[disc:]

fftsq = np.zeros((maxim,len(t_new)))
fftsq_c = np.zeros((maxim,len(t_new)))

Omegas = list()
    
for i in range(maxim):
    Omega = Omega + 0.5
    Omegas.append(Omega)

    ivp_sol = solve_ivp(system_ivp,[0,t0], X0, method=method, rtol=1e-6, t_eval=t, dense_output=False, events=None, vectorized=True)
    sol1 = ivp_sol.y
    t1 = ivp_sol.t
    sol = np.zeros([t_len,4])
    for I in [0,1,2,3]:
         sol[:,I] = np.interp(t, t1, sol1[I,:])
    
    solx = sol[:,0]
    solx = solx[disc:] # discarding the transient period
    soldx = sol[:,1]
    soly = sol[:,2]
    soly = soly[disc:] # discarding the transient period
    soldy = sol[:,3]      
    
    # FFT 
    fftx = np.fft.fft(solx)
    fftx2 = (fftx * np.conjugate(fftx))**0.5
    ffty = np.fft.fft(soly)
    ffty2 =(ffty * np.conjugate(ffty))**0.5

    fftsq[i,:] = (ffty2**2 + fftx2**2)**0.5
    
    fftsq[i,:] = np.log(fftsq_c[i,:])   


# In[19]:


from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

# auxiliary for correct scale for plotting
end_point = max([Omega,omega_nat])+1
x_scale = np.arange(0,end_point,2*pi/(t_new[-1]-t_new[0]))

T = len(x_scale)
end = 25

fig = plt.figure(figsize=(30,20))
ax = Axes3D(fig)

ax.view_init(elev=80, azim=-90)

X = (np.ones((T,end)).T*x_scale)
Y = (np.zeros((T,end))+ omega).T
Z = fftsq[:end,:T]

# Plot a basic wireframe.
ax.set_title('Waterfall plot of FFTs')
ax.plot_wireframe(X,Y,Z, rstride=1, cstride=0)
ax.grid('off')
ax.set_ylabel(ylabel='Driving frequency'+r' $\Omega$'+'\n \n ')
ax.set_xlabel(xlabel='\n \n Frequency')
ax.plot3D(xs=[omega_nat,omega_nat],ys=[1.15,6.1],zs=[8,8],lw=5,c='C1',label='Natural frequency')
ax.plot3D(xs=[1.15,1.15+100*0.05],ys=[1.15,6.1],zs=[8,8],lw=5,c='C2',label='Driving frequency')

#ax.set_axis_off()
#ax.set_frame_on(True)
ax.set_zticks(ticks=[])
plt.legend()
# plt.savefig('Im2/waterfall_VDH_biggap_with_lines.png')
plt.show()


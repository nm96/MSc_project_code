# A file for recording parameter value sets that lead to interesting results.

# 'Transitional chaos' in the VdH model:
eps = 0.18 # Rotor eccentricity
Om = 3.1 # Driving frequency
m = 10 # Mass (per unit length)
c = 0.05 # Damping coefficient
k = 1 # Stiffness coefficient
h = 1 # Gap width
k_c = 50 # Stator stiffness parameter for VdH model

# Parameter values from Flores 2009
h = 0.2*10**-3 # Clearance (m)
b = 40*10**-3 # Bearing length (m)
R2 = 9.8*10**-3 # Journal Radius (m)
mu = 4 # Oil Viscosity (Pascals) - at 40C 
m = 0.13 # Journal mass (kg) - big ?? on this one..
k = 10 # ??
c = 1 # ??
eps = 0.3
Om = 4.1

beta = 12*mu*b*R2**3/h**2

print("Flores Beta = {:.4f}".format(beta))


# Rich spectrum of subharmonics in the NHS lubrication model:
eps = 0.2 # Rotor eccentricity
Om = 4.1 # Driving frequency
m = 10 # Mass (per unit length)
c = 0.5 # Damping coefficient
k = 10 # Stiffness coefficient
h = 1 # Gap width
mu = 10**-7 # Viscosity
b = 1 # Bearing length
R2 = 100 # Radius

beta = 12*mu*b*R2**3/h**2

print("Critical Beta = {:.4f}".format(beta))

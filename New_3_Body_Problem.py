## Imports
import astropy.constants as constants
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as ode
## Constants
# Gravitational Constant
G = constants.G.value

# Masses
m_earth = constants.M_earth.value
m_sun = constants.M_sun.value
m_moon = 7.34767309e22

# Distances in astronomic unit
r_sun = np.array([0,0,0])
r_earth = np.array([1, 0, 0])
r_moon = np.array([1.002569, 0, 0])

# Instantaneous speed of bodies (in units of 2pi/365/24/60/60) - angular momentum
v_sun = np.array([0,0,0])
v_earth = np.array([0, 1, 0])
v_moon = np.array([0, 0.0747945, 0]) #27.3/365

## Variables
m1, m2, m3 = 1, 1, 1
r1, r2, r3 = [0,0,0], [1,1,1], [2,2,2]

## Math
def coupled_de(r, t, G, m1, m2, m3):
    # separate r & v
    r1 = r[:3]
    r2 = r[3:6]
    r3 = r[6:9]
    v1 = r[9:12]
    v2 = r[12:15]
    v3 = r[15:]

    # Differences
    r12 = np.linalg.norm(r1-r2)
    r23 = np.linalg.norm(r2-r3)
    r13 = np.linalg.norm(r1-r3)

    # second derivatives
    a1 = -m2*(r1-r2)/np.abs(r12)**3 - m3*(r1-r3)/np.abs(r13)**3
    a2 = -m3*(r2-r3)/np.abs(r23)**3 + m1*(r1-r2)/np.abs(r12)**3
    a3 = +m1*(r1-r3)/np.abs(r13)**3 + m2*(r2-r3)/np.abs(r23)**3
    # first derivatives
    u1 = v1
    u2 = v2
    u3 = v3
    args = (a1, a2, a3, u1, u2, u3)
    all = np.concatenate(args)
    return all

params = np.array([r_sun, r_earth, r_moon, v_sun, v_earth, v_moon]) #Initial parameters
params = params.flatten() #Flatten to make 1D array
t = np.linspace(0,20,500)
three_body_sol = ode.odeint(coupled_de, params, t, args=(G,m_sun,m_earth,m_moon))

r1_sol=three_body_sol[:,:3]
r2_sol=three_body_sol[:,3:6]
r3_sol=three_body_sol[:,6:9]
## Graphing
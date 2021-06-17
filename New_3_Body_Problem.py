## Imports
import astropy.constants as constants
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as ode
## Constants
# Gravitational Constant
G = constants.G

# Masses
m_earth = constants.M_earth.value
m_sun = constants.M_sun
m_moon = 7.34767309e22

# Distances in astronomic unit
r_sun = np.array([0,0,0])
r_earth = np.array([1, 0, 0])
r_moon = np.array([1+0.002569, 0, 0])

# Instantaneous speed of bodies (in units of 2pi/365/24/60/60) - angular momentum
v_sun = np.array([0,0,0])
v_earth = np.array([0, 1, 0])
v_moon = np.array([0, 27.3/365, 0])

## Variables
m1, m2, m3 = 1, 1, 1
r1, r2, r3 = [0,0,0], [1,1,1], [2,2,2]

## Math
def coupled_de(m, r1, r2, r3, G, t):

    # Differences
    r12 = np.linalg.norm(r1-r2)
    r23 = np.linalg.norm(r2-r3)
    r13 = np.linalg.norm(r1-r3)

    # second derivatives
    ddr1dr = -G*m[1]*r12/np.abs(r12)**3 - G*m[2]*r13/np.abs(r13)**3
    ddr2dr = -G*m[2]*r23/np.abs(r23)**3 + G*m[0]*r12/np.abs(r12)**3
    ddr3dr = +G*m[0]*r13/np.abs(r13)**3 + G*m[1]*r23/np.abs(r23)**3

    args = (ddr1dr, ddr2dr, ddr3dr)
    all = np.concatenate(args)
    return all

params = np.array([r_sun, r_earth, r_moon, v_sun, v_earth, v_moon]) #Initial parameters
params = params.flatten() #Flatten to make 1D array
t = np.linspace(0,20,500)
three_body_sol = ode.odeint(coupled_de,params,t,args=(G,m_sun,m_earth,m_moon))

r1_sol=three_body_sol[:,:3]
r2_sol=three_body_sol[:,3:6]
r3_sol=three_body_sol[:,6:9]
## Graphing
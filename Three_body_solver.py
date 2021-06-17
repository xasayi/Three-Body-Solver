import matplotlib.pyplot as plt
import numpy as np
import scipy as sci
import scipy.integrate as ode
#from celluloid import Camera
import mpl_toolkits
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation


#constants
G = 6.67408e-11#gravitational constant
m = 1.989e+30 #kg #mass of the sun
r = 148e9 #m #distance between sun and moon
v = 30000 #m/s #relative velocity of earth around the sun
t = 365*24*3600 #s #orbital period of earth
K1=G*t*m/(r**2*v)
K2=v*t/r

flat = 0
round = 0
other1 = 0
other2 = 1
other3 = 0

if round==1:
    m1 = 1                   #sun
    m2 = 5.972e24/1.989e30   #earth
    m3 = 0.07346e24/1.989e30 #moon
    r1 = np.array([0, 0, 0])      #initial position of sun
    r2 = np.array([-1, 0, 0])     #initial position of earth
    r3 = np.array([-1.001, 0, 0]) #initial position of moon
    v1 = np.array([0, 0, 0])
    v2 = np.array([0, 0, -1])
    v3 = np.array([0.05, 00, -1.05])

if flat==1:
    m1 = 1                   #sun
    m2 = 5.972e24/1.989e30   #earth
    m3 = 0.07346e24/1.989e30 #moon
    r1 = np.array([0, 0.00000321868, 0.00000321868])     #initial position of sun
    r2 = np.array([0, 0, 0])      #initial position of earth
    r3 = np.array([0, -0.00000321868, 0.00000321868])#initial position of moon
    v1 = np.array([0, 0, -1])
    v2 = np.array([0, 0, 0])
    v3 = np.array([0.05, 0, -1.05])
if other1==1:
    m1 = 1                   #sun
    m2 = 1   #earth
    m3 = 0.001 #moon
    r1 = np.array([0, 0, 0])     #initial position of sun
    r2 = np.array([-1, 0, 0])      #initial position of earth
    r3 = np.array([-2, 0, 0])#initial position of moon
    v1 = np.array([0, 0, -1])
    v2 = np.array([0, 0, 0])
    v3 = np.array([0, 0, -1])
if other2==1:
    m1 = 1                   #sun
    m2 = 0.001   #earth
    m3 = 1 #moon
    r1 = np.array([0, 0, 0])     #initial position of sun
    r2 = np.array([-1, 0, 0])      #initial position of earth
    r3 = np.array([-1.5, 0, 0])#initial position of moon
    v1 = np.array([0, 0, -1])
    v2 = np.array([1, 0, 0])
    v3 = np.array([0, 0, -1])
if other3==1:
    m1 = 1                   #sun
    m2 = 1   #earth
    m3 = 1 #moon
    r1 = np.array([0, 0, 0])     #initial position of sun
    r2 = np.array([-2, 0, 0])      #initial position of earth
    r3 = np.array([-4, 0, 0])#initial position of moon
    v1 = np.array([0, -1, 0])
    v2 = np.array([0.05, 0, 0])
    v3 = np.array([0, 0, -1])

def ThreeBodyEquations(w,t,G,m1,m2,m3):
    r1 = w[:3]
    r2 = w[3:6]
    r3 = w[6:9]
    v1 = w[9:12]
    v2 = w[12:15]
    v3 = w[15:18]
    r12 = np.linalg.norm(r2-r1) #vector from earth to sun
    r13 = np.linalg.norm(r3-r1) #vector from earth to moon
    r23 = np.linalg.norm(r3-r2) #vector from moon to sun

    dv1ydt = K1*m2*(r2-r1)/r12**3 + K2*m3*(r3-r1)/r13**3
    dv2ydt = K1*m1*(r1-r2)/r12**3 + K2*m3*(r3-r2)/r23**3
    dv3ydt = K1*m1*(r1-r3)/r13**3 + K2*m2*(r2-r3)/r23**3
    dr1ydt = K2*v1
    dr2ydt = K2*v2
    dr3ydt = K2*v3
    args = (dr1ydt, dr2ydt, dr3ydt, dv1ydt, dv2ydt, dv3ydt)
    der = np.concatenate(args)
    return der

init_params = np.array([r1,r2,r3,v1,v2,v3]) #Initial parameters
init_params = init_params.flatten() #Flatten to make 1D array
time_span = np.linspace(0,20,500)
three_body_sol = ode.odeint(ThreeBodyEquations,init_params,time_span,args=(G,m1,m2,m3))

r1_sol=three_body_sol[:,:3]
r2_sol=three_body_sol[:,3:6]
r3_sol=three_body_sol[:,6:9]

##
#Create figure
fig = plt.figure()
ax = fig.add_subplot(111,projection="3d")
#camera = Camera(fig)

def update_lines(num, dataLines, lines):
    for line, data in zip(lines, dataLines):
        line.set_data(data[0:2, :num])
        line.set_3d_properties(data[2, :num])
    return lines

#fig = plt.figure()
ax = p3.Axes3D(fig)

# Lines to plot in 3D
data = np.array([[r1_sol[:,0],r1_sol[:,1],r1_sol[:,2]],[r3_sol[:,0],r3_sol[:,1],r3_sol[:,2]],[r2_sol[:,0],r2_sol[:,1],r2_sol[:,2]]])
color = ['gold','darkgrey','skyblue']
label = ['sun', 'moon', 'earth']
s = [1, 3, 1]
lines = [ax.plot(data[i][0, 0:1], data[i][1, 0:1], data[i][2, 0:1], color[i], label = label[i], linewidth = s[i] )[0] for i in [0,1,2]]
#lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]

ax.scatter(r1_sol[-1,0],r1_sol[-1,1],r1_sol[-1,2],color="gold",marker="*",s=100)

ax.set_xlim(-5,5)
ax.set_ylim(-5,5)
ax.set_zlim(-5,5)

ax.set_title("Sun Moon Earth System\n")
ax.legend(loc="upper left")

line_ani = animation.FuncAnimation(fig, update_lines, 50, fargs=(data, lines),
                                   interval=50, blit=True, repeat=True)
plt.show()
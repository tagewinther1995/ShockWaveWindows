import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from pyctp import ljs_bh
import os
import imageio
from mpl_toolkits.mplot3d import axes3d
# Modify system path
import sys
sys.path.append('../pycThermopack/')
# Importing pyThermopack
ljs = ljs_bh.ljs_bh()
ljs.init()
ljs.set_tmin(temp=2.0)

# Get parameters for Argon
sigma, eps = ljs.get_sigma_eps()
# For konversjon
m = 6.6335209e-26  # kg
# Plot phase envelope using a1-a3
z = np.array([1])

N = 100  # number of points to discretize
L = 1.0
xc = 0.25
v = 1
X = np.linspace(0, L, N) # position along the rod
h = L / (N - 1)
radius = 0.2 #[M]
hVolume = np.pi*radius*radius*h
molarMass = 0.039948




def odefunc(u, t):
    dudt = np.zeros(2*len(X))
    u[N] = 0
    u[N-1] = 1

    dudt[0] =  - np.exp(5.73*(u[0]-u[N]))+np.exp(-11.46*(u[0]-u[N]))
    dudt[-1] = np.exp(5.73*(u[N-1]-u[-1]))-np.exp(-11.46*(u[N-1]-u[-1]))

    #Central differencing at points on lfs and rhs
    #u1 at 1 and N-2
    #w = 12
    #dudx1 = (-u[3]+16*u[2]-30*u[1]+16*u[0])/(w*h*h)
    dudx1 = (u[2]-2*u[1]+u[0])/(h*h)
    dudt[1] = 0.024* dudx1 - np.exp(5.73*(u[1]-u[N+1]))+np.exp(-11.46*(u[1]-u[N+1]))

    #dudx1 = (16*u[N-1]-30*u[N-2]+16*u[N-3]-u[N-4])/(w*h*h)
    dudx2 = (u[N-1]-2*u[N-2]+u[N-3])/(h*h)
    dudt[N-2] = 0.024* dudx1 - np.exp(5.73*(u[N-2]-u[-2]))+np.exp(-11.46*(u[N-2]-u[-2]))

    #dudx2 = (-u[N+1+2]+16*u[N+1+1]-30*u[N+1]+16*u[N])/(w*h*h)
    dudx2 = (u[N+2]-2*u[N+1]+u[N])/(h*h)
    dudt[N+1] = 0.170* dudx2 + np.exp(5.73*(u[1]-u[N+1]))+np.exp(-11.46*(u[1]-u[N+1]))

    #dudx2 = (16*u[-1]-30*u[-2]+16*u[-3]-u[-4])/(w*h*h)
    dudx2 = (u[-3]-2*u[-2]+u[-1])/(h*h)
    dudt[-2] = 0.170* dudx2 + np.exp(5.73*(u[N-2]-u[-2]))+np.exp(-11.46*(u[N-2]-u[-2]))

    #forward and backward differencing at points on lfs and rhs
    #forward
    #dudx1 = (2*u[1]-5*u[2]+4*u[3]-u[4])/(h*h)
    #dudt[1] = 0.024* dudx1 - np.exp(5.73*(u[1]-u[N+1]))+np.exp(-11.46*(u[1]-u[N+1]))

    #dudx2 = (2*u[N+1]-5*u[N+2]+4*u[N+3]-u[N+4])/(h*h)
    #dudt[N+1] = 0.170* dudx2 + np.exp(5.73*(u[1]-u[N+1]))+np.exp(-11.46*(u[1]-u[N+1]))

    #backward

    #dudx1 = (2*u[N-2]-5*u[N-3]+4*u[N-4]-u[N-5])/(h*h)
    #dudt[N-2] = 0.024* dudx1 - np.exp(5.73*(u[N-2]-u[-2]))+np.exp(-11.46*(u[N-2]-u[-2]))

    #dudx2 = (2*u[-2]-5*u[-3]+4*u[-4]-u[-5])/(h*h)
    #dudt[-2] = 0.170* dudx2 + np.exp(5.73*(u[N-2]-u[-2]))+np.exp(-11.46*(u[N-2]-u[-2]))

    # now for the internal nodes
    for i in range(2, N-2):
        

        #du1dx = (u[i+1]-2*u[i]+u[i-1])/(h*h)
        #du2dx = (u[N+i+1]-2*u[N+i]+u[N+i-1])/(h*h)
        #higher order
        du1dx = (-u[i+2]+16*u[i+1]-30*u[i]+16*u[i-1]-u[i-2])/(12*h*h)
        du2dx = (-u[N+i+2]+16*u[N+i+1]-30*u[N+i]+16*u[N+i-1]-u[N+i-2])/(12*h*h)

        #dudt[i] = 0.024* (u[i+1]-2*u[i]+u[i-1])/(h*h) - np.exp(5.73*(u[i]-u[N+i]))+np.exp(-11.46*(u[i]-u[N+i]))
        dudt[i] = 0.024* du1dx - np.exp(5.73*(u[i]-u[N+i]))+np.exp(-11.46*(u[i]-u[N+i]))
        
        #dudt[i+N] = 0.170* (u[N+i+1]-2*u[N+i]+u[N+i-1])/(h*h) + np.exp(5.73*(u[i]-u[N+i]))-np.exp(-11.46*(u[i]-u[N+i]))

        dudt[i+N] = 0.170* du2dx + np.exp(5.73*(u[i]-u[N+i]))-np.exp(-11.46*(u[i]-u[N+i]))
        
    
    
    return dudt

init = []
realSol =[]
for i in range(len(X)):
    u0 = 1
    init.append(u0)
for i in range(len(X)):
    init.append(0)


tspan = np.linspace(0.0, 2, 100)
sol = odeint(odefunc, init, tspan)
filenames = []
u1 = []
u2 = []
lastPoint = sol[-1][:N]
for i in range(0, len(tspan)):
    plt.plot(X, sol[i][N:], 'g.',label= "u2")
    plt.plot(X,sol[i][:N],'bx', label = "u1")
    u2.append(sol[i][N:])
    u1.append(sol[i][:N])
    plt.legend()

    #plt.show()

    #stringName = "Python/figures/modelling/" + str(i) + 'timeStep.png'
    #filenames.append(stringName)
    #plt.savefig(stringName)
    plt.clf()

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, projection='3d')
u1 = np.array(u1)
u2 = np.array(u2)
# Plot a 3D surface


X,tspan = np.meshgrid(X,tspan)
ax.plot_surface(X, tspan, u1)
plt.show()
plt.savefig("Python/figures/modelling/3Du2plot.PNG")







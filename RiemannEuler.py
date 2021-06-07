import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from pyctp import ljs_bh
import os
import imageio
from scipy import constants
# Modify system path
import sys

def evenlySpaced(start,stop,N):
    liste = [start]
    a = start
    for i in range(N):
        b = a+ (stop / (N - 1))
        liste.append(b)
        a = a +(stop/(N-1))

    return liste
sys.path.append('../pycThermopack/')
# Importing pyThermopack
ljs = ljs_bh.ljs_bh()
ljs.init()
ljs.set_tmin(temp=2.0)

NA = 6.02214076e23
# Get parameters for Argon
sigma, eps = ljs.get_sigma_eps()
aeps = eps*constants.k
# For konversjon
m = 6.6335209e-26  # kg
# Plot phase envelope using a1-a3

z = np.array([1])
# Time
tStop = 0.25 #s
K = 40
tspan = np.linspace(0.0, tStop, K)

# Space
N = 750# number of points to discretize
L = 1.0
Lratio = 16

X = np.linspace(-L, L, N)  # position along the rod

T = 124  # K


h = L / (N - 1)
k = tStop/(K-1)
radius = L/Lratio  # [m]
hVolume = np.pi*radius*radius*h  # [m3]
molarMass = 0.039948  # kg/mol

rho0 = 1
rhoL = 2
vL = 0
rhoR = 1
vR = 0
c = 1
dxdt = h/k
gamma = c*c/rho0



def exactSolution(xList, tList, rho0, rhoL, vL, rhoR, vR, c):
    rhoMatrix = []
    vMatrix = []
    #The Common denominator
    rhoC = 2*c*rho0

    a1 =  (c*rhoL - rho0*vL)/rhoC
    a2 =  (c*rhoL + rho0*vL)/rhoC

    b1 = (c*rhoR - rho0*vR)/rhoC
    b2 = (c*rhoR + rho0*vR)/rhoC


    #initial
    rhoMatrix.append([])
    vMatrix.append([])
    for j in range(len(xList)):
        if xList[j] <= 0:
            rhoMatrix[0].append(rhoL)
            vMatrix[0].append(vL)
        else:
            rhoMatrix[0].append(rhoR)
            vMatrix[0].append(vR)


    for i in range(1,len(tList)):
        rhoMatrix.append([])
        vMatrix.append([])
        for j in range(len(xList)):
            if tList[i] <= (-xList[j]/c) and tList[i] > 0:
                rho = a1*rho0 + a2*rho0
                rhoMatrix[i].append(rho)
                v = -a1*c + a2*c
                vMatrix[i].append(v)
            elif (abs(xList[j]/c)) >= 0 and (abs(xList[j]/c))< tList[i]:
                rho = b1*rho0 + a2*rho0
                rhoMatrix[i].append(rho)
                v = -b1*c + a2*c
                vMatrix[i].append(v)
            elif tList[i] > 0 and tList[i] <= (xList[j]/c):
                rho = b1*rho0 + b2*rho0
                rhoMatrix[i].append(rho)
                v = -b1*c + b2*c
                vMatrix[i].append(v)
    return [rhoMatrix, vMatrix]


def plottingFunctionHeat(yArraySecond, yArrayFourth, yLabel, name, animate=True, plot2order=False):
    filenames = []
    for i in range(0, len(yArraySecond)):
        
        if plot2order == False:
            plt.plot(X, yArrayFourth[i], "c-", label="MOL: 4.th order central differencing")
        plt.plot(X, yArraySecond[i], color="darkred",
                     linestyle="solid", label="Analytical Solution")
        plt.legend()
        #plt.plot(X, sol[i][N:])
        plt.ylabel(yLabel)
        plt.xlabel("x [m]")
        stringName = "Python/figures/modelling/" + str(i) + 'timeStep.png'
        filenames.append(stringName)
        plt.savefig(stringName)
        plt.clf()

    if animate == True:
        gifName = "Python/figures/modelling/" + name + "EulerMovie.gif"
        with imageio.get_writer(gifName, mode='I', duration=1) as writer:
            for filename in filenames:
                image = imageio.imread(filename)
                writer.append_data(image)

        for m in filenames:
            os.remove(m)

    return None


def plottingFunctionHeat4figuresDensity(yArraySecond, yArrayFourth, yArrayForce,yArrayLax, yLabel, name, animate=True, plot2order=False):
    filenames = []
    for i in range(0, len(yArraySecond)):
        plt.rcParams['axes.spines.right'] = False
        plt.rcParams['axes.spines.top'] = False
        plt.axis([X[0], X[-1], 0, (max(yArraySecond[0])+max(yArraySecond[0])*0.4)])
        if plot2order == False:
            plt.plot(X, yArrayFourth[i], color = "c", linestyle="dotted", label="MOL: 4.th order central differencing")
            plt.plot(X,yArrayLax[i],color = 'orange',linestyle="dashed", label = "Lax-Friedrich")
            plt.plot(X, yArrayForce[i], color ="darkgreen",linestyle="dashdot",  label="FORCE")
        plt.plot(X, yArraySecond[i], color="darkred",linestyle="solid", label="Analytical")

        plt.legend(loc = 1)
        #plt.plot(X, sol[i][N:])
        plt.ylabel(yLabel)
        plt.xlabel("x [m]")
        stringName = "Python/figures/modelling/" + str(tspan[i]) + 'timeStep.png'
        filenames.append(stringName)
        plt.savefig(stringName)
        plt.clf()

    if animate == True:
        gifName = "Python/figures/modelling/" + name + "EulerMovie.gif"
        with imageio.get_writer(gifName, mode='I', duration=1) as writer:
            for filename in filenames:
                image = imageio.imread(filename)
                writer.append_data(image)

        for m in filenames:
            os.remove(m)

    return None

def plottingFunctionHeat4figuresVelocity(yArraySecond, yArrayFourth, yArrayForce,yArrayLax, yLabel, name, animate=True, plot2order=False):
    filenames = []
    for i in range(0, len(yArraySecond)):
        plt.rcParams['axes.spines.right'] = False
        plt.rcParams['axes.spines.top'] = False
        plt.axis([X[0], X[-1], 0, (max(yArraySecond[1])+max(yArraySecond[1])*0.8)])
        if plot2order == False:
            plt.plot(X, yArrayFourth[i], color = "c", linestyle="dotted", label="MOL: 4.th order central differencing")
            plt.plot(X,yArrayLax[i],color = 'orange',linestyle="dashed", label = "Lax-Friedrich")
            plt.plot(X, yArrayForce[i], color ="darkgreen",linestyle="dashdot",  label="FORCE")
        plt.plot(X, yArraySecond[i], color="darkred",linestyle="solid", label="Analytical")

        plt.legend(loc = 1)
        #plt.plot(X, sol[i][N:])
        plt.ylabel(yLabel)
        plt.xlabel("x [m]")
        stringName = "Python/figures/modelling/" + str(tspan[i]) + 'timeStep.png'
        filenames.append(stringName)
        plt.savefig(stringName)
        plt.clf()

    if animate == True:
        gifName = "Python/figures/modelling/" + name + "EulerMovie.gif"
        with imageio.get_writer(gifName, mode='I', duration=1) as writer:
            for filename in filenames:
                image = imageio.imread(filename)
                writer.append_data(image)

        for m in filenames:
            os.remove(m)

    return None


def odefunc1(u, t):

    dudt = np.zeros(2*len(X))

    dudt[0] = 0
    dudt[N-1] = 0

    dudt[N] = 0
    dudt[-1] = 0
    #u[N] = 0
    #u[-1] =  0

    #dndt[-1] = 0
    #dvdt[-1] = 0
    #dudt[-1] = (u[-1] - u[-2])/h

    # now for the internal nodes
    for i in range(1, N-1):

        #dudt[i] = - (u[i + 1]*u[N+i + 1] - u[i - 1]*u[N+i -1]) / (2*h)
        dudt[i] = -rho0 *(u[N+i + 1] - u[N+i - 1]) / (2*h)
        dudt[N+i] = -(c*c/rho0)*((u[i + 1] - u[i - 1]) / (2*h))

    return dudt


def odefunc1LaxFORCE(u, t):

    dudt = np.zeros(2*len(X))

    dudt[0] = 0
    dudt[N-1] = 0

    dudt[N] = 0
    dudt[-1] = 0
    #u[N] = 0
    #u[-1] =  0

    #dndt[-1] = 0
    #dvdt[-1] = 0
    #dudt[-1] = (u[-1] - u[-2])/h

    # now for the internal nodes
    for i in range(1, N-1):

        #dudt[i] = - (u[i + 1]*u[N+i + 1] - u[i - 1]*u[N+i -1]) / (2*h)
        #dudt[i] = -rho0 *(u[N+i + 1] - u[N+i - 1]) / (2*h)
        #dudt[N+i] = -(c*c/rho0)*((u[i + 1] - u[i - 1]) / (2*h))
        FnPlus = 0.5*rho0*(u[N+i]+ u[N+i+1])+ 0.5*dxdt*(u[i]-u[i+1])
        FvPlus = 0.5*gamma*(u[i]+ u[i+1])+ 0.5*dxdt*(u[N+i]-u[N+i+1])

        FnMinus = 0.5*rho0*(u[N+i-1]+ u[N+i])- 0.5*dxdt*(u[i]-u[i-1])
        FvMinus = 0.5*gamma*(u[i-1]+ u[i])- 0.5*dxdt*(u[N+i]-u[N+i-1])
    
        dudt[i] = -(FnPlus-FnMinus)/h
        dudt[N+i] = -(FvPlus-FvMinus)/h


    return dudt


def odefunc1FORCE(u, t):

    dudt = np.zeros(2*len(X))

    dudt[0] = 0
    dudt[N-1] = 0

    dudt[N] = 0
    dudt[-1] = 0
    #u[N] = 0
    #u[-1] =  0

    #dndt[-1] = 0
    #dvdt[-1] = 0
    #dudt[-1] = (u[-1] - u[-2])/h

    # now for the internal nodes
    for i in range(1, N-1):

        #dudt[i] = - (u[i + 1]*u[N+i + 1] - u[i - 1]*u[N+i -1]) / (2*h)
        #dudt[i] = -rho0 *(u[N+i + 1] - u[N+i - 1]) / (2*h)
        #dudt[N+i] = -(c*c/rho0)*((u[i + 1] - u[i - 1]) / (2*h))

        #i + 1/2
        FnLfPlus = 0.5*rho0*(u[N+i]+ u[N+i+1])+ 0.5*dxdt*(u[i]-u[i+1])
        FvLfPlus = 0.5*gamma*(u[i]+ u[i+1])+ 0.5*dxdt*(u[N+i]-u[N+i+1])
        QnRiPlus = 0.5*(u[i+1]+u[i])+0.5*rho0*(1/dxdt)*(u[N+i]-u[N+i+1])
        QvRiPlus = 0.5*(u[N+i+1]+u[N+i])+0.5*gamma*(1/dxdt)*(u[i]-u[i+1])
        FnRiPlus = rho0*QvRiPlus
        FvRiPlus = gamma*QnRiPlus
        FnPlus = 0.5*(FnLfPlus+FnRiPlus)
        FvPlus = 0.5*(FvLfPlus+FvRiPlus)

        #i -1/2

        FnLf = 0.5*rho0*(u[N+i-1]+ u[N+i])- 0.5*dxdt*(u[i]-u[i-1])
        FvLf = 0.5*gamma*(u[i-1]+ u[i])- 0.5*dxdt*(u[N+i]-u[N+i-1])
        #inncorrect fill inn
        QnRi = 0.5*(u[i-1]+u[i])-0.5*rho0*(1/dxdt)*(u[N+i]-u[N+i-1])
        QvRi = 0.5*(u[N+i-1]+u[N+i])-0.5*gamma*(1/dxdt)*(u[i]-u[i-1])

        #QnRi = 0.5*(u[i-1]+u[i])-0.5*rho0*(1/dxdt)*(u[N+i]-u[N+i-1])
        #QvRi = 0.5*(u[N+i-1]+u[N+i])-0.5*gamma*(1/dxdt)*(u[i]-u[i-1])

        FnRi = rho0*QvRi
        FvRi = gamma*QnRi
        FnMinus = 0.5*(FnLf+FnRi)
        FvMinus = 0.5*(FvLf+FvRi)
        dudt[i] = -(FnPlus-FnMinus)/h
        dudt[N+i] = -(FvPlus-FvMinus)/h


    return dudt


def odefunc1spatial(u, x):

    dudx = np.zeros(2*len(X))

    dudx[0] = 0
    dudx[K-1] = 0

    dudx[K] = 0
    dudx[-1] = 0
    #u[N] = 0
    #u[-1] =  0

    #dndt[-1] = 0
    #dvdt[-1] = 0
    #dudt[-1] = (u[-1] - u[-2])/h

    # now for the internal nodes
    for i in range(1, K-1):

        #dudt[i] = - (u[i + 1]*u[N+i + 1] - u[i - 1]*u[N+i -1]) / (2*h)
        dudx[i] = -1/rho0 *(u[K+i+1] - u[K+i]) / (k)
        dudx[K+i] = -(rho0/(c*c))*((u[i + 1] - u[i]) / (k))

    return dudx


def odefunc1Alt(u, t):

    dudt = np.zeros(2*len(X))

    dudt[0] = 0
    dudt[N-1] = 0

    dudt[N] = 0
    dudt[-1] = 0
    #u[N] = 0
    #u[-1] =  0
    # discretization at point 1 and N-2
    #dudt[i] = - (u[i + 1]*u[N+i + 1] - u[i - 1]*u[N+i -1]) / (2*h)
    # discretization at point 1( rho[1], rhov[N+1])
    dudt[1] = - rho0*(u[N+1 + 1] - u[N+1-1]) / (2*h)
    # discretization at point N-2, rho[N-2] rhov [-2]
    dudt[N-2] = - rho0*(u[-1] - u[-3]) / (2*h)

    dudt[N+1] = -(c*c/rho0)*((u[2] - u[0]) / (2*h))
    
    dudt[-2] = -(c*c/rho0)*((u[-1] - u[-3]) / (2*h))

    # now for the internal nodes
    for i in range(2, N-2):

        #dudt[i] = - (u[i + 1]*u[N+i + 1] - u[i - 1]*u[N+i -1]) / (2*h)
        #dudt[i] = - (u[N+i + 1] - u[N+i -1]) / (2*h)
        dudt[i] = - rho0*(-u[N+i + 2]+8*u[N+i + 1] - 8 *u[N+i - 1] + u[N+i - 2]) / (12*h)
        dudt[i+N] = -(c*c/rho0)*(-u[i + 2]+8*u[i + 1] -
                  8*u[i - 1] + (u[i - 2])) / (12*h)

    return dudt

def isotermMatrixMaker(solutions):
    densities = []
    velocities = []
    for i in range(len(solutions)):
        densities.append([])
        velocities.append([])
        for j in range(N):
            densities[i].append(solutions[i][j])
            velocity = solutions[i][N+j]
            velocities[i].append(velocity)
    return [densities, velocities]

def isotermMatrixMakerSpatial(solutions):
    densities = []
    velocities = []
    for i in range(len(solutions)):
        densities.append([])
        velocities.append([])
        for j in range(K):
            densities[i].append(solutions[i][j])
            velocity = solutions[i][N+j]
            velocities[i].append(velocity)
    return [densities, velocities]

#Initial conditions
# rho0, rhoL, vL, rhoR, vR, c

init1 = []
init2 = []
for j in range(len(X)):
    if X[j] <= 0:
        init1.append(rhoL)
        init2.append(vL)
    else:
        init1.append(rhoR)
        init2.append(vR)
init  = init1 + init2

exactSol = exactSolution(X,tspan, rho0,rhoL,vL,rhoR,vR,c)
sol1 = odeint(odefunc1, init, tspan)
sol2 = odeint(odefunc1FORCE, init, tspan)
sol3 = odeint(odefunc1LaxFORCE, init, tspan)
#sol1 = odeint(odefunc1spatial,init, X)

numSol1 = isotermMatrixMaker(sol1)
numSol2 = isotermMatrixMaker(sol2)
numSol3 = isotermMatrixMaker(sol3)

numDensity1 = numSol1[0]
numVelocity1 = numSol1[1]
numDensity2 = numSol2[0]
numVelocity2 = numSol2[1]
numDensity3 = numSol3[0]
numVelocity3 = numSol3[1]

density = exactSol[0]
velocity = exactSol[1]
plottingFunctionHeat4figuresDensity(density, numDensity1,numDensity2,numDensity3, "ρ [kg/m3]", "excactSolDensity",True, False)
plottingFunctionHeat4figuresVelocity(velocity, numVelocity1, numVelocity2, numVelocity3,"v [m/s]", "excactSolVelocity",True, False)

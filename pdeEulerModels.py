import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from pyctp import ljs_bh
import os
import imageio
from scipy import constants
# Modify system path
import sys
sys.path.append('../pycThermopack/')
# Importing pyThermopack
ljs = ljs_bh.ljs_bh()
ljs.init()
ljs.set_tmin(temp=2.0)

NA = 6.02214076e23
m = 6.6335209e-26


def calc_reduced_T(Ta, eps):
    """ Calculate reduced temperature
    """
    Tstar = np.zeros_like(Ta)
    Tstar = Ta/eps
    return Tstar

def calc_reduced_time(tReal, eps,sigma):
    """ Calculate reduced temperature
    """
    Tstar = tReal*(1/sigma)*np.sqrt(eps/m)
    
    return Tstar


def calc_real_T(Tstar, eps):

    # Calcualte real temperature
    Ta = np.zeros_like(Tstar)
    Ta = Tstar * eps
    return Ta


def calc_reduced_rho(rhoa, sigma):
    """ Calculate reduced density
    """
    rhoStar = np.zeros_like(rhoa)
    rhoStar = sigma**3*NA*rhoa
    return rhoStar


def calc_reduced_sos(velocity, m, eps):
    # Calculate reduced speed of sound
    sosStar = np.zeros_like(velocity)
    sosStar = velocity*np.sqrt(m/eps)
    return sosStar


def calc_Real_P(pStar, sigma, eps):

    pR = np.zeros_like(pStar)
    pR = pStar * eps/(sigma**3)

    return pR


def calc_real_h(hStar, m, eps):
    hR = np.zeros_like(hStar)
    hR = hStar * eps/(m)
    return hR

def calc_real_u(uStar, m, eps):
    uR = 0.039948*uStar * eps/(m)
    return uR


def calc_reduced_h(hr, m, eps):
    #hStar = np.zeros_like(hr[0])
    h = hr/0.039948  # Change from mole to mass basis
    hStar = h * m/eps
    return hStar


def calc_reduced_s(sr, m):

    s = sr[0]/0.039948
    sStar = s*m/constants.k
    return sStar


def calc_reduced_u(Ur, m, eps, Np):
    U = Ur[0]/(m*Np)
    Ustar = U*m/eps
    return Ustar


def calc_Real_volume(Vstar, sigma):
    Vr = Vstar*sigma**3
    return Vr


def calc_reduced_p(pR, sigma, eps):
    #pStar = (pR[0]/eps)*sigma**3
    pStar = (pR/eps)*sigma**3
    return pStar


def listToArray(yList, rows, columns):
    npList = np.array(yList)
    shape = (rows, columns)
    npList.reshape(shape)
    return npList






# Get parameters for Argon
sigma, eps = ljs.get_sigma_eps()
aeps = eps*constants.k
# For konversjon
  # kg
# Plot phase envelope using a1-a3

z = np.array([1])
#tStop = 0.00005
#ikke gå under 25
K = 25

# Space
N = 25# number of points to discretize
nemdScale = True

if nemdScale == True:
    #Nemd case
    Lreduced = 150#*sigma
    L = 150*sigma #[m]
    #tStop = 3e-12
    tStop = 3e-15
    #xc = L/2
    #L = 1.0
    Xreduced = np.linspace(0, Lreduced, N)
else:
    L = 1
    Xreduced = np.linspace(0, L, N)
    tStop = 0.00005



T = 124  # K
X = np.linspace(0, L, N)
tspan = np.linspace(0.0, tStop, K)

h = L / (N - 1)
k = tStop/(K-1)
Lratio = 32
radius = L/Lratio  # [m]

hVolume = np.pi*radius*radius*h  # [m3]
#hVolume = h*radius*radius
molarMass = 0.039948  # kg/mol
dxdt = h/k


def normalDistriubtionMaker(xValue,centre,width):
    return((1/(width*np.sqrt(2*np.pi))*np.exp(-0.5*((xValue-centre)/width)**2)))


def calc_real_rho_kg(rhoStar, sigma):
    rhoa = rhoStar/(sigma**3*NA)  # mol/m3
    rhoa = rhoa*molarMass  # mol/n3 * kg/mol
    return rhoa


def calc_reduced_rho_kg(rhoa, sigma):
    """ Calculate reduced density
    """
    rhoStar = rhoa*(sigma**3*NA)  # kg/m3*m3*mol-1
    rhoStar = rhoStar/molarMass  # kg/mol / kg/mol
    return rhoStar


def calc_reduced_sos(velocity, m, eps):
    # Calculate reduced speed of sound
    sosStar = velocity*np.sqrt(m/eps)
    return sosStar

def calc_real_u_kg(uStar, m, eps):
    uR = uStar * eps/(m)
    return uR


def pressureIsotermCalc(densityMatrix, temperature, realPressure=True):
    pressureMatrix = []
    for i in range(int(K)):
        pressureMatrix.append([])
        for j in range(N):
            if realPressure == True:
                moles = (hVolume*densityMatrix[i][j])/molarMass
                p = ljs.pressure_tv(temperature, hVolume, [moles])[0]
                pressureMatrix[i].append(p)
            else:
                moles = (hVolume*densityMatrix[i][j])/molarMass
                p = ljs.pressure_tv(temperature, hVolume, [moles])[0]
                p = calc_reduced_p(p, sigma, aeps)
                pressureMatrix[i].append(p)
    return pressureMatrix


def viscosityIsotermCalc(densityMatrix, Temperature):
    visMatrix = []
    for i in range(int(K)):
        visMatrix.append([])
        for j in range(N):
            visco = argonViscosity(Temperature, densityMatrix[i][j])
            visMatrix[i].append(visco)
    return visMatrix


def isotermMatrixMaker(solutions):
    densities = []
    velocities = []
    for i in range(len(solutions)):
        densities.append([])
        velocities.append([])
        for j in range(N):
            densities[i].append(solutions[i][j])
            velocity = solutions[i][N+j]/solutions[i][j]
            velocities[i].append(velocity)
    return [densities, velocities]


def energyMatrixMaker(solutions):
    densities = []
    velocities = []
    internalEnergies = []
    for i in range(len(solutions)):
        densities.append([])
        velocities.append([])
        internalEnergies.append([])
        for j in range(N):
            densities[i].append(solutions[i][j])
            velocity = solutions[i][N+j]/solutions[i][j]
            velocities[i].append(velocity)
            internalEnergy = solutions[i][2*N+j]/solutions[i][j] - 0.5*solutions[i][N+j]**2/(solutions[i][j]**2)
            internalEnergies[i].append(internalEnergy)
    return [densities, velocities, internalEnergies]


def isotermMatrixMakerDimless(solutions):
    densities = []
    velocities = []
    for i in range(len(solutions)):
        densities.append([])
        velocities.append([])
        for j in range(N):
            rhoStar = calc_reduced_rho_kg(solutions[i][j], sigma)
            densities[i].append(rhoStar)
            velocity = solutions[i][N+j]/solutions[i][j]
            vStar = calc_reduced_sos(velocity, m, aeps)
            velocities[i].append(vStar)
    return [densities, velocities]


def argonViscosity(temp, rho):
    #rho [kg/m3]
    # convert density to correct dimensions
    rhor = rho/(molarMass*1000)  # mol/dm3
    # dilute gas viscosity
    ek = 143.2
    sig = 0.335
    Tc = 150.687  # kelvin
    rhoc = 13.40743  # mol/dm3
    b = [0.431, -0.4623, 0.08406, 0.005341, -0.00331]
    Tdot = temp/ek
    aList = []
    omega = 1
    for i in range(len(b)):
        a = np.exp(b[i]*(np.log(Tdot)**(i)))
        aList.append(a)
    for i in range(len(b)):
        omega = omega*aList[i]
    etaDilute = 0.0266958*np.sqrt(temp*39.948)/(omega*sig**2)  # mu * Pa*s
    # Residual viscosity
    Ni = [12.19, 13.99, 0.005027, -18.93, -6.698, -3.827]
    ti = [0.42, 0, 0.95, 0.5, 0.9, 0.8]
    di = [1, 2, 10, 5, 1, 2]
    li = [0, 0, 0, 2, 4, 4]
    fi = [0, 0, 0, 1, 1, 1]
    tau = Tc/temp
    delta = rhor/rhoc
    etr = 0
    for i in range(len(Ni)):
        eti = Ni[i]*tau**(ti[i])*delta**(di[i])*np.exp(-fi[i]*delta**(li[i]))
        etr = etr + eti

    eta = (etaDilute + etr) * 10**(-6)  # Pa * s
    return eta


def plottingFunction(yArraySecond, yArrayFourth, yLabel, name, animate=True, plot2order=False, dimless = False):
    filenames = []
    for i in range(0, len(yArraySecond)):
        if min(yArraySecond[1]) >= 0:
            #plt.axis([0, X[-1], 0, (max(yArraySecond[0])+max(yArraySecond[0])*0.4)])
            #plt.axis([0, X[-1], 0, (max(yArraySecond[0])+max(yArraySecond[0])*0.4)])
            print("ech")
        else:
            #plt.axis([0, X[-1], min(yArraySecond[i]) + 0.25*min(yArraySecond[i]), (max(yArraySecond[i])+max(yArraySecond[i])*0.25)])
            #plt.axis([0, X[-1], min(yArraySecond[1]) + 0.25*min(yArraySecond[1]),(max(yArraySecond[1])+max(yArraySecond[1])*0.25)])
            print("ech")

        if dimless == False:
            if plot2order == False:
                plt.plot(X, yArrayFourth[i], color="darkred",
                        linestyle="dashed", label="MOL: Lax-Friedrich")
            plt.plot(X, yArraySecond[i], color="dodgerblue",
                    linestyle="solid", label="MOL: FORCE")
            plt.xlabel("x [m]")
        else:
            if plot2order == False:
                plt.plot(Xreduced, yArrayFourth[i], color="darkred",
                        linestyle="dashed", label="MOL: 4.th order central differencing")
            plt.plot(Xreduced, yArraySecond[i], color="dodgerblue",
                    linestyle="solid", label="MOL: FORCE")
            plt.xlabel("x*")
        plt.legend(loc = 2)
        #plt.plot(X, sol[i][N:])
        plt.ylabel(yLabel)
        
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



def plottingFunctionHeat(yArraySecond, yArrayFourth, yLabel, name, animate=True, plot2order=False):
    filenames = []
    for i in range(0, len(yArraySecond)):
        if min(yArraySecond[i]) > 0:
            #plt.axis([0, X[-1], 0, (max(yArraySecond[i])+max(yArraySecond[i])*0.25)])
            #plt.axis([0, X[-1], 0, (max(yArraySecond[1])+max(yArraySecond[1])*0.25)])
            print("ech")
        else:
            #plt.axis([0, X[-1], min(yArraySecond[i]) + 0.25*min(yArraySecond[i]), (max(yArraySecond[i])+max(yArraySecond[i])*0.25)])
            #plt.axis([0, X[-1], min(yArraySecond[1]) + 0.25*min(yArraySecond[1]),(max(yArraySecond[1])+max(yArraySecond[1])*0.25)])
            print("ech")

        plt.plot(X, yArraySecond[i], "c--", label="second order")
        if plot2order == False:
            plt.plot(X, yArrayFourth[i], color="darkred",
                     linestyle="solid", label="Fourth order")
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

def plotting3figuresDensity(yCentre, yLax, yForce,yLabel, name, animate=True, dimless = False):
    filenames = []
    snapShotsMoments = [0,3,5,10,15,20]
    for i in range(0, len(yCentre)):
        
        if dimless == True:
            plt.rcParams['axes.spines.right'] = False
            plt.rcParams['axes.spines.top'] = False
            plt.axis([Xreduced[0], Xreduced[-1], 0, (max(yCentre[0])+max(yCentre[0])*0.4)])
            plt.plot(Xreduced, yCentre[i], color = "c", linestyle="dotted", label="MOL:4.th order central differencing")
            plt.plot(Xreduced,yLax[i],color = 'orange',linestyle="dashed", label = "Lax-Friedrich")
            plt.plot(Xreduced, yForce[i], color ="darkgreen",linestyle="dashdot",  label="FORCE")
            plt.legend(loc = 1)
            plt.ylabel(yLabel)
            plt.xlabel("x*")
        else:
            plt.rcParams['axes.spines.right'] = False
            plt.rcParams['axes.spines.top'] = False
            plt.axis([X[0], X[-1], 0, (max(yCentre[0])+max(yCentre[0])*0.4)])
            plt.plot(X, yCentre[i], color = "c", linestyle="dotted", label="MOL:4.th order central differencing")
            plt.plot(X,yLax[i],color = 'orange',linestyle="dashed", label = "Lax-Friedrich")
            plt.plot(X, yForce[i], color ="darkgreen",linestyle="dashdot",  label="FORCE")
            plt.legend(loc = 1)
            plt.ylabel(yLabel)
            plt.xlabel("x*")
        stringName = "Python/figures/modelling/" + str(tspan[i]) + 'timeStep.png'
        filenames.append(stringName)
        if i in snapShotsMoments:
            shotname = "Python/figures/modelling/snapshots/" + name + str(calc_reduced_time(tspan[i],aeps,sigma)) +'.png'
            plt.savefig(shotname)
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

def plotting3figuresVelocity(yCentre, yLax, yForce,yLabel, name, animate=True, dimless = False):
    filenames = []
    snapShotsMoments = [0,3,5,10,15,20]
    for i in range(0, len(yCentre)):
        
        if dimless == True:
            plt.rcParams['axes.spines.right'] = False
            plt.rcParams['axes.spines.top'] = False
            plt.axis([Xreduced[0], Xreduced[-1], (min(yCentre[4])+min(yCentre[4])*0.4), (max(yCentre[4])+max(yCentre[4])*0.4)])
            plt.plot(Xreduced, yCentre[i], color = "c", linestyle="dotted", label="MOL:4.th order central differencing")
            plt.plot(Xreduced,yLax[i],color = 'orange',linestyle="dashed", label = "Lax-Friedrich")
            plt.plot(Xreduced, yForce[i], color ="darkgreen",linestyle="dashdot",  label="FORCE")
            plt.legend(loc = 2)
            plt.ylabel(yLabel)
            plt.xlabel("x*")
        else:
            plt.rcParams['axes.spines.right'] = False
            plt.rcParams['axes.spines.top'] = False
            plt.axis([X[0], X[-1], (min(yCentre[4])+min(yCentre[4])*0.4), (max(yCentre[4])+max(yCentre[4])*0.4)])
            plt.plot(X, yCentre[i], color = "c", linestyle="dotted", label="MOL:4.th order central differencing")
            plt.plot(X,yLax[i],color = 'orange',linestyle="dashed", label = "Lax-Friedrich")
            plt.plot(X, yForce[i], color ="darkgreen",linestyle="dashdot",  label="FORCE")
            plt.legend(loc = 2)
            plt.ylabel(yLabel)
            plt.xlabel("x*")
        
        stringName = "Python/figures/modelling/" + str(tspan[i]) + 'timeStep.png'
        filenames.append(stringName)
        if i in snapShotsMoments:
            shotname = "Python/figures/modelling/snapshots/" + name + str(tspan[i]) +'.png'
            plt.savefig(shotname)
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


def odefunc1Visco(u, t):

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
        dudt[i] = - (u[N+i + 1] - u[N+i - 1]) / (2*h)
        moles1 = (hVolume*u[i+1])/molarMass  # m3*kg/m3 / (kg/mol)
        moles0 = (hVolume*u[i-1])/molarMass
        dPdx = (ljs.pressure_tv(T, hVolume, [moles1])[
                0] - (ljs.pressure_tv(T, hVolume, [moles0])[0])) / (2*h)
        #dnnvdx = (u[i + 1]*u[N+i + 1]**2 - u[i - 1]*u[N+i -1]**2) / (2*h)
        dnnvdx = (u[N+i + 1]**(2)/u[i + 1] -
                  (u[N+i - 1]**(2)/u[i - 1])) / (2*h)
        if i == 1:
            # setter V(i-2) == 0
            dviscodx = (argonViscosity(
                T, u[i+1])*(u[N+i + 2]/u[i + 2] - u[N+i]/u[i]) - argonViscosity(T, u[i-1])*(u[N+i]/u[i]))/(3*h*h)
        elif i == N-2:
            # setter V(i+2) == 0
            dviscodx = (argonViscosity(T, u[i+1])*(u[N+i]/u[i]) - argonViscosity(
                T, u[i-1])*(u[N+i]/u[i] - u[N+i-2]/u[i-2]))/(3*h*h)
        else:
            dviscodx = (argonViscosity(T, u[i+1])*(u[N+i + 2]/u[i + 2] - u[N+i]/u[i]) -
                        argonViscosity(T, u[i-1])*(u[N+i]/u[i] - u[N+i-2]/u[i-2]))/(3*h*h)

        dudt[i+N] = -(dPdx + dnnvdx) + dviscodx

    return dudt


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
        dudt[i] = - (u[N+i + 1] - u[N+i - 1]) / (2*h)
        moles1 = (hVolume*u[i+1])/molarMass  # m3*kg/m3 / (kg/mol)
        moles0 = (hVolume*u[i-1])/molarMass
        dPdx = (ljs.pressure_tv(T, hVolume, [moles1])[
                0] - (ljs.pressure_tv(T, hVolume, [moles0])[0])) / (2*h)
        #dnnvdx = (u[i + 1]*u[N+i + 1]**2 - u[i - 1]*u[N+i -1]**2) / (2*h)
        dnnvdx = (u[N+i + 1]**(2)/u[i + 1] -
                  (u[N+i - 1]**(2)/u[i - 1])) / (2*h)
    

        dudt[i+N] = -(dPdx + dnnvdx)

    return dudt

def odefunc1LaxFORCE(u, t):

    dudt = np.zeros(2*len(X))
    #Once again we're lazy and assume the values are constant at the boundaries
    dudt[0] = 0
    dudt[N-1] = 0

    dudt[N] = 0
    dudt[-1] = 0
    
    for i in range(1, N-1):
        molesR = (hVolume*u[i+1])/molarMass #i+1
        moles0 = (hVolume*u[i])/molarMass #i
        molesL = (hVolume*u[i-1])/molarMass #i-1
        pR = ljs.pressure_tv(T, hVolume, [molesR])[0]
        p0 = ljs.pressure_tv(T, hVolume, [moles0])[0]
        pL = ljs.pressure_tv(T, hVolume, [molesL])[0]
        FnPlus = 0.5*(u[N+i]+ u[N+i+1])+ 0.5*dxdt*(u[i]-u[i+1])
        FnvPlus = 0.5*(u[N+i]*u[N+i]/u[i] +p0 +u[N+i+1]*u[N+i+1]/u[i+1]+ pR)+ 0.5*dxdt*(u[N+i]-u[N+i+1])
        #Lax-Friedrich scheme
        FnLf = 0.5*(u[N+i-1]+ u[N+i])- 0.5*dxdt*(u[i]-u[i-1])
        FnvLf = 0.5*(u[N+i-1]*u[N+i-1]/u[i-1] +pL +u[N+i]*u[N+i]/u[i]+ p0)- 0.5*dxdt*(u[N+i]-u[N+i-1])
        #inncorrect fill inn

        #QnRi = 0.5*(u[i-1]+u[i])-0.5*(1/dxdt)*(u[N+i]-u[N+i-1])
        #QnvRi = 0.5*(u[N+i-1]+u[N+i])-0.5*(1/dxdt)*(-u[N+i-1]*u[N+i-1]/u[i-1] -pL +u[N+i]*u[N+i]/u[i]+ p0)

        #FnRi = QnvRi
        #FnvRi = QnvRi*QnvRi/QnRi + 0.5*(p0+pL)
        #FnMinus = 0.5*(FnLf+FnRi)
        #FnvMinus = 0.5*(FnvLf+FnvRi)
        FnMinus = FnLf
        FnvMinus = FnvLf
        dudt[i] = -(FnPlus-FnMinus)/h
        dudt[N+i] = -(FnvPlus-FnvMinus)/h

    return dudt

def odefunc1FORCE(u, t):

    dudt = np.zeros(2*len(X))
    #Once again we're lazy and assume the values are constant at the boundaries
    dudt[0] = 0
    dudt[N-1] = 0

    dudt[N] = 0
    dudt[-1] = 0
    
    for i in range(1, N-1):
        molesR = (hVolume*u[i+1])/molarMass #i+1
        moles0 = (hVolume*u[i])/molarMass #i
        molesL = (hVolume*u[i-1])/molarMass #i-1
        pR = ljs.pressure_tv(T, hVolume, [molesR])[0]
        p0 = ljs.pressure_tv(T, hVolume, [moles0])[0]
        pL = ljs.pressure_tv(T, hVolume, [molesL])[0]
        #i + 0.5

        FnLfPlus = 0.5*(u[N+i]+ u[N+i+1])+ 0.5*dxdt*(u[i]-u[i+1])
        FnvLfPlus = 0.5*(u[N+i]*u[N+i]/u[i] +p0 +u[N+i+1]*u[N+i+1]/u[i+1]+ pR) + 0.5*dxdt*(u[N+i]-u[N+i+1])
        QnRiPlus = 0.5*(u[i+1]+u[i])+0.5*(1/dxdt)*(u[N+i]-u[N+i+1])
        QnvRiPlus = 0.5*(u[N+i+1]+u[N+i])+0.5*(1/dxdt)*(-u[N+i+1]*u[N+i+1]/u[i+1] -pR +u[N+i]*u[N+i]/u[i]+ p0)
        FnRiPlus = QnvRiPlus
        FnvRiPlus = QnvRiPlus*QnvRiPlus/QnRiPlus + 0.5*(p0+pR)
        FnPlus = 0.5*(FnLfPlus+FnRiPlus)
        FnvPlus = 0.5*(FnvLfPlus+FnvRiPlus)

        #i -1/2
        FnLf = 0.5*(u[N+i-1]+ u[N+i])- 0.5*dxdt*(u[i]-u[i-1])
        FnvLf = 0.5*(u[N+i-1]*u[N+i-1]/u[i-1] +pL +u[N+i]*u[N+i]/u[i]+ p0)- 0.5*dxdt*(u[N+i]-u[N+i-1])
        QnRi = 0.5*(u[i-1]+u[i])-0.5*(1/dxdt)*(u[N+i]-u[N+i-1])
        QnvRi = 0.5*(u[N+i-1]+u[N+i])-0.5*(1/dxdt)*(-u[N+i-1]*u[N+i-1]/u[i-1] -pL +u[N+i]*u[N+i]/u[i]+ p0)

        FnRi = QnvRi
        FnvRi = QnvRi*QnvRi/QnRi + 0.5*(p0+pL)
        FnMinus = 0.5*(FnLf+FnRi)
        FnvMinus = 0.5*(FnvLf+FnvRi)
        dudt[i] = -(FnPlus-FnMinus)/h
        dudt[N+i] = -(FnvPlus-FnvMinus)/h

    return dudt   


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
    dudt[1] = - (u[N+1 + 1] - u[N+1-1]) / (2*h)
    moles1 = (hVolume*u[1+1])/molarMass
    moles0 = (hVolume*u[1-1])/molarMass
    dPdx = (ljs.pressure_tv(T, hVolume, [moles1])[
            0] - (ljs.pressure_tv(T, hVolume, [moles0])[0])) / (2*h)
    #dnnvdx = (u[i + 1]*u[N+i + 1]**2 - u[i - 1]*u[N+i -1]**2) / (2*h)
    dnnvdx = (u[N+1 + 1]**(2)/u[1 + 1] - (u[N+1 - 1]**(2)/u[1 - 1])) / (2*h)
    dudt[1+N] = -(dPdx + dnnvdx)
    # discretization at point N-2, rho[N-2] rhov [-2]
    dudt[N-2] = - (u[-1] - u[-3]) / (2*h)
    moles1 = (hVolume*u[N-1])/molarMass
    moles0 = (hVolume*u[N-3])/molarMass
    dPdx = (ljs.pressure_tv(T, hVolume, [moles1])[0] - (ljs.pressure_tv(T, hVolume, [moles0])[0])) / (2*h)
    dnnvdx = (u[-1]**(2)/u[N - 1] - u[-3]**(2)/u[N - 3]) / (2*h)
    dudt[-2] = -(dPdx + dnnvdx)

    # now for the internal nodes
    for i in range(2, N-2):

        #dudt[i] = - (u[i + 1]*u[N+i + 1] - u[i - 1]*u[N+i -1]) / (2*h)
        #dudt[i] = - (u[N+i + 1] - u[N+i -1]) / (2*h)
        dudt[i] = - (-u[N+i + 2]+8*u[N+i + 1] - 8 *
                     u[N+i - 1] + u[N+i - 2]) / (12*h)
        moles2 = (hVolume*u[i+2])/molarMass
        moles1 = (hVolume*u[i+1])/molarMass
        molesN1 = (hVolume*u[i-1])/molarMass
        molesN2 = (hVolume*u[i-2])/molarMass
        dPdx = (-ljs.pressure_tv(T, hVolume, [moles2])[0] + 8*ljs.pressure_tv(T, hVolume, [moles1])[0] - 8*(
            ljs.pressure_tv(T, hVolume, [molesN1])[0]) + ljs.pressure_tv(T, hVolume, [molesN2])[0]) / (12*h)
        #dnnvdx = (u[i + 1]*u[N+i + 1]**2 - u[i - 1]*u[N+i -1]**2) / (2*h)
        dnnvdx = (-u[N+i + 2]**(2)/u[i + 2]+8*u[N+i + 1]**(2)/u[i + 1] -
                  8*(u[N+i - 1]**(2)/u[i - 1]) + (u[N+i - 2]**(2)/u[i - 2])) / (12*h)


                
        dudt[i+N] = -(dPdx + dnnvdx)

    return dudt


def odefunc2(u, t):
    #Energy equation
    dudt = np.zeros(3*len(X))
    # Boundary conditions of rho
    dudt[0] = 0
    dudt[N-1] = 0
    # Boundary condition of rho*v
    dudt[N] = 0
    dudt[2*N-1] = 0
    #Heating at end of tube case:
    # Boundary condition of E = rho*(u)
    dudt[2*N] = 0
    dudt[-1] = 0
    #Heating at end of tube case:
    # now for the internal nodes
    for i in range(1, N-1):
        #u[i] = rho[i]
        #u[N+i] = rho*v
        #u[2N+i] = E = rho(u + 0.5*v^2)

        #dudt[i] = - (u[i + 1]*u[N+i + 1] - u[i - 1]*u[N+i -1]) / (2*h)
        # Determining Temperature and pressure at necessary points
        #
        internalEnergy1 = molarMass*(u[2*N+i+1]/u[i+1] - 0.5*u[N+i+1]**(2)/(u[i+1]**2)) #J/mol
        internalEnergy0 = molarMass*(u[2*N+i-1]/u[i-1] - 0.5*u[N+i-1]**(2)/(u[i-1]**2)) #J/mol
        specificVolume1 = molarMass/u[i+1] # kg/mol *M3/kg
        specificVolume0 = molarMass/u[i-1] # kg/mol *M3/kg
        tempPress1 = ljs.two_phase_uvflash(z, internalEnergy1, specificVolume1)
        press1 = tempPress1[1]#Pa
        tempPress0 = ljs.two_phase_uvflash(z, internalEnergy0, specificVolume0)
        press0 = tempPress0[1] #pa
        # Continuity equation
        dudt[i] = - (u[N+i + 1] - u[N+i - 1]) / (2*h)
        # derivatives for momentum equation
        dPdx = (press1 - press0) / (2*h)
        #dnnvdx = (u[i + 1]*u[N+i + 1]**2 - u[i - 1]*u[N+i -1]**2) / (2*h)
        dnnvdx = (u[N+i + 1]**(2)/u[i + 1] -(u[N+i - 1]**(2)/u[i - 1])) / (2*h)
        # The momentum equation
        dudt[i+N] = -(dPdx + dnnvdx)
        #Derivatives for the heat equation
        #dEvdx = E[i+1]*v[i+1] - E[i+1]*v[i+1]
        dEvdx = ((u[2*N+i+1]*u[N+i + 1]/u[i + 1]) -(u[2*N+i-1]*u[N+i - 1]/u[i - 1]))/(2*h)
        dPvdx = ((press1*u[N+i + 1]/u[i + 1]) -(press0*u[N+i - 1]/u[i - 1]))/(2*h)
        dudt[2*N+i] = -(dEvdx + dPvdx)
    return dudt


def odefunc2FORCE(u, t):
    #Energy equation
    dudt = np.zeros(3*len(X))
    # Boundary conditions of rho
    dudt[0] = 0
    dudt[N-1] = 0
    # Boundary condition of rho*v
    dudt[N] = 0
    dudt[2*N-1] = 0
    #Heating at end of tube case:
    # Boundary condition of E = rho*(u)
    dudt[2*N] = 0
    dudt[-1] = 0
    #Heating at end of tube case:
    # now for the internal nodes
    for i in range(1, N-1):
        #u[i] = rho[i]
        #u[N+i] = rho*v
        #u[2N+i] = E = rho(u + 0.5*v^2)

        #dudt[i] = - (u[i + 1]*u[N+i + 1] - u[i - 1]*u[N+i -1]) / (2*h)
        # Determining Temperature and pressure at necessary points
        #
        internalEnergyR = molarMass*(u[2*N+i+1]/u[i+1] - 0.5*u[N+i+1]**(2)/(u[i+1]**2)) #J/mol
        internalEnergy0 = molarMass*(u[2*N+i]/u[i] - 0.5*u[N+i]**(2)/(u[i]**2)) #J/mol
        internalEnergyL = molarMass*(u[2*N+i-1]/u[i-1] - 0.5*u[N+i-1]**(2)/(u[i-1]**2)) #J/mol
        specificVolumeR = molarMass/u[i+1] # kg/mol *M3/kg
        specificVolume0 = molarMass/u[i]
        specificVolumeL = molarMass/u[i-1] # kg/mol *M3/kg
        tempPress1 = ljs.two_phase_uvflash(z, internalEnergyR, specificVolumeR)
        pR = tempPress1[1]#Pa
        tempPress0 = ljs.two_phase_uvflash(z, internalEnergy0, specificVolume0)
        p0= tempPress0[1] #pa
        tempPressL = ljs.two_phase_uvflash(z, internalEnergyL, specificVolumeL)
        pL = tempPressL[1] #pa
        eR = (u[2*N+i+1]+pR)*u[N+i + 1]/u[i + 1]
        e0 = (u[2*N+i]+p0)*u[N+i]/u[i]
        eL = (u[2*N+i-1]+pL)*u[N+i - 1]/u[i -1]
        
        FnLfPlus = 0.5*(u[N+i]+ u[N+i+1])+ 0.5*dxdt*(u[i]-u[i+1])
        FnvLfPlus = 0.5*(u[N+i]*u[N+i]/u[i] +p0 +u[N+i+1]*u[N+i+1]/u[i+1]+ pR)+ 0.5*dxdt*(u[N+i]-u[N+i+1])
        FeLfPlus = 0.5*(e0 +eR)+ 0.5*dxdt*(u[2*N+i]-u[2*N+i+1])

        QnRiPlus = 0.5*(u[i+1]+u[i])+0.5*(1/dxdt)*(u[N+i]-u[N+i+1])
        QnvRiPlus = 0.5*(u[N+i+1]+u[N+i])+0.5*(1/dxdt)*(-u[N+i+1]*u[N+i+1]/u[i+1] -pR +u[N+i]*u[N+i]/u[i]+ p0)
        QeRiPLus = 0.5*(u[2*N +i+1]+u[2*N + i])+0.5*(1/dxdt)*(e0-eL)

        FnRiPlus = QnvRiPlus
        FnvRiPlus = QnvRiPlus*QnvRiPlus/QnRiPlus + 0.5*(p0+pR)
        FeRiPlus = (QeRiPLus +0.5*(p0+pR))*QnvRiPlus/QnRiPlus

        FnPlus = 0.5*(FnLfPlus+FnRiPlus)
        FnvPlus = 0.5*(FnvLfPlus+FnvRiPlus)
        FePlus = 0.5*(FeLfPlus+ FeRiPlus)
        #Lax-Friedrich scheme
        FnLf = 0.5*(u[N+i-1]+ u[N+i])- 0.5*dxdt*(u[i]-u[i-1])
        FnvLf = 0.5*(u[N+i-1]*u[N+i-1]/u[i-1] +pL +u[N+i]*u[N+i]/u[i]+ p0)- 0.5*dxdt*(u[N+i]-u[N+i-1])
        FeLf = 0.5*(eL + e0)- 0.5*dxdt*(u[2*N+i]-u[2*N+i-1])
            

        QnRi = 0.5*(u[i-1]+u[i])-0.5*(1/dxdt)*(u[N+i]-u[N+i-1])
        QnvRi = 0.5*(u[N+i-1]+u[N+i])-0.5*(1/dxdt)*(-u[N+i-1]*u[N+i-1]/u[i-1] -pL +u[N+i]*u[N+i]/u[i]+ p0)
        QeRi = 0.5*(u[2*N +i-1]+u[2*N + i])-0.5*(1/dxdt)*(e0-eL)

        FnRi = QnvRi
        FnvRi = QnvRi*QnvRi/QnRi + 0.5*(p0+pL)
        FeRi = (QeRi +0.5*(p0+pL))*QnvRi/QnRi

        FnMinus = 0.5*(FnLf+FnRi)
        FnvMinus = 0.5*(FnvLf+FnvRi)
        FeMinus = 0.5*(FeLf+ FeRi)

        dudt[i] = -(FnPlus-FnMinus)/h
        dudt[N+i] = -(FnvPlus-FnvMinus)/h
        dudt[2*N+i] = -(FePlus-FeMinus)/h
    return dudt



def odefunc2FORCELAX(u, t):
    #Energy equation
    dudt = np.zeros(3*len(X))
    # Boundary conditions of rho
    dudt[0] = 0
    dudt[N-1] = 0
    # Boundary condition of rho*v
    dudt[N] = 0
    dudt[2*N-1] = 0
    #Heating at end of tube case:
    # Boundary condition of E = rho*(u)
    dudt[2*N] = 0
    dudt[-1] = 0
    #Heating at end of tube case:
    # now for the internal nodes
    for i in range(1, N-1):
        #u[i] = rho[i]
        #u[N+i] = rho*v
        #u[2N+i] = E = rho(u + 0.5*v^2)

        #dudt[i] = - (u[i + 1]*u[N+i + 1] - u[i - 1]*u[N+i -1]) / (2*h)
        # Determining Temperature and pressure at necessary points
        #
        internalEnergyR = molarMass*(u[2*N+i+1]/u[i+1] - 0.5*u[N+i+1]**(2)/(u[i+1]**2)) #J/mol
        internalEnergy0 = molarMass*(u[2*N+i]/u[i] - 0.5*u[N+i]**(2)/(u[i]**2)) #J/mol
        internalEnergyL = molarMass*(u[2*N+i-1]/u[i-1] - 0.5*u[N+i-1]**(2)/(u[i-1]**2)) #J/mol
        specificVolumeR = molarMass/u[i+1] # kg/mol *M3/kg
        specificVolume0 = molarMass/u[i]
        specificVolumeL = molarMass/u[i-1] # kg/mol *M3/kg
        tempPress1 = ljs.two_phase_uvflash(z, internalEnergyR, specificVolumeR)
        pR = tempPress1[1]#Pa
        tempPress0 = ljs.two_phase_uvflash(z, internalEnergy0, specificVolume0)
        p0= tempPress0[1] #pa
        tempPressL = ljs.two_phase_uvflash(z, internalEnergyL, specificVolumeL)
        pL = tempPressL[1] #pa
        eR = (u[2*N+i+1]+pR)*u[N+i + 1]/u[i + 1]
        e0 = (u[2*N+i]+p0)*u[N+i]/u[i]
        eL = (u[2*N+i-1]+pL)*u[N+i - 1]/u[i -1]
        
        FnPlus = 0.5*(u[N+i]+ u[N+i+1])+ 0.5*dxdt*(u[i]-u[i+1])
        FnvPlus = 0.5*(u[N+i]*u[N+i]/u[i] +p0 +u[N+i+1]*u[N+i+1]/u[i+1]+ pR)+ 0.5*dxdt*(u[N+i]-u[N+i+1])
        FePlus = 0.5*(e0 +eR)+ 0.5*dxdt*(u[2*N+i]-u[2*N+i+1])

    
        #Lax-Friedrich scheme
        FnMinus = 0.5*(u[N+i-1]+ u[N+i])- 0.5*dxdt*(u[i]-u[i-1])
        FnvMinus = 0.5*(u[N+i-1]*u[N+i-1]/u[i-1] +pL +u[N+i]*u[N+i]/u[i]+ p0)- 0.5*dxdt*(u[N+i]-u[N+i-1])
        FeMinus = 0.5*(eL + e0)- 0.5*dxdt*(u[2*N+i]-u[2*N+i-1])
            

        dudt[i] = -(FnPlus-FnMinus)/h
        dudt[N+i] = -(FnvPlus-FnvMinus)/h
        dudt[2*N+i] = -(FePlus-FeMinus)/h
    return dudt


def main():
    #Testing area for uvflash
    uReal = calc_real_u_kg(-1.493908,m,aeps)
    dens = calc_real_rho_kg(0.6, sigma)
    moles = (hVolume*dens)/molarMass
    uThermo = ljs.internal_energy_tv(124,hVolume,[moles])
    uThermo = uThermo[0]/moles
    
    molarDensity = dens/molarMass
    molarVolume = 1/molarDensity

   
    tempTest = ljs.two_phase_uvflash(z, uThermo, molarVolume)
    phase = tempTest[-1]
    a = ljs.get_phase_type(phase)
    pressureList = []
    init = []
    realSol = []
    isoterm = False
    test = ljs.two_phase_uvflash(z, 1334.2859740744148, 6.746923211402392e-05)
    #make them unitless
    reducedU = calc_reduced_h(1333.2859740744148, m, aeps)
    density = molarMass/6.746923211402392e-05
    reducedRho = calc_reduced_rho_kg(density,sigma)


    a = test[0]
    print(a)
    print(test[1])
    # print(test)
    if isoterm == True:
        for i in range(len(X)):
            #if X[i] >= ((L*0.5)-(L*0.02)) and X[i] <= ((L*0.5)+(L*0.02)):


                #rho = calc_real_rho_kg(0.9, sigma)
                # rho = 0.0006 #[kg/m3]
                #init.append(rho)
            #else:
                #rho = calc_real_rho_kg(0.6, sigma)
                #rho = 0.0004
                #init.append(rho)
            if nemdScale == True:
                rho = np.exp(-0.005*(Xreduced[i]-Lreduced/2)**2) +0.6
            else:
                rho = np.exp(-20*(Xreduced[i]-Lreduced/2)**2) +0.6
            #rho = normalDistriubtionMaker(X[i], L/2,L/10) + 0.6
            rho = calc_real_rho_kg(rho, sigma)
            init.append(rho)
        plt.plot(X,init)
        plt.show()
        
        for i in range(len(X)):
            init.append(0*init[i])

        sol1 = odeint(odefunc1Alt, init, tspan)
        sol2 = odeint(odefunc1LaxFORCE, init, tspan)
        sol3 = odeint(odefunc1FORCE,init,tspan)
        
        

        solutions = isotermMatrixMaker(sol1)
        densities1 = solutions[0]
        velocities1 = solutions[1]

        #viscosity = viscosityIsotermCalc(densities1, T)

        solutions = isotermMatrixMaker(sol2)
        densities2 = solutions[0]
        velocities2 = solutions[1]

        solutions = isotermMatrixMaker(sol3)
        densities3 = solutions[0]
        velocities3 = solutions[1]

        filenames = []
        pressures1 = pressureIsotermCalc(densities1, T)
        pressures2 = pressureIsotermCalc(densities2, T)
        pressures3 = pressureIsotermCalc(densities3, T)
        pressures1Dim = pressureIsotermCalc(densities1, T, False)
        pressures2Dim = pressureIsotermCalc(densities2, T, False)
        pressures3Dim = pressureIsotermCalc(densities3, T, False)
            
        # real units
        #plottingFunction(densities1, densities2, "n [kg/m3]", "force100density", True)
        #plottingFunction(velocities1, velocities2, "v [m/s]", "force100velocity", True)
        #plottingFunction(pressures1, pressures2, "P [pa]", "force100pressure", True)
        #plottingFunction(viscosity, pressures1,"mu [Pa*s]", "500viscosity", True, True)
        plotting3figuresDensity(densities1, densities2, densities3, "n [kg/m3]", "TRUEforce1000density", True)
        plotting3figuresVelocity(velocities1, velocities2, velocities3, "v [m/s]", "TRUEforce1000velocity", True)
        plotting3figuresDensity(pressures1, pressures2, pressures3 ,"P [pa]", "TRUEforce1000pressure", True)
        # Dimensionless units

        solutions = isotermMatrixMakerDimless(sol1)
        densities1 = solutions[0]
        velocities1 = solutions[1]

        solutions = isotermMatrixMakerDimless(sol2)
        densities2 = solutions[0]
        velocities2 = solutions[1]

        solutions = isotermMatrixMakerDimless(sol3)
        densities3 = solutions[0]
        velocities3 = solutions[1]

        plotting3figuresDensity(densities1, densities2, densities3, "n*", "TRUEforce1000Dimlessdensity", True,True)
        plotting3figuresVelocity(velocities1, velocities2, velocities3, "v*", "TRUEforce1000Dimlessvelocity", True, True)
        plotting3figuresDensity(pressures1Dim, pressures2Dim, pressures3Dim ,"P*", "TRUEforce1000Dimlesspressure", True, True)


        #plottingFunction(densities1, densities2, "n*", "force100DimlessDensity", True,False,True)
        #plottingFunction(velocities1, velocities2, "v*", "force100DimlessVelocity", True,False,True)
        #plottingFunction(pressures1Dim, pressures2Dim,"P*", "force100fDimlessPressure", True,False,True)
    else:
        temp = []
        dens = calc_real_rho_kg(0.6, sigma)
        moles = (hVolume*dens)/molarMass
        v0 = 0 #[m/s]
        #initializing density, density*velocity and total energy  E = density*U + 0.5*density*velocity**2
        for i in range(len(X)):
            rho = calc_real_rho_kg(0.6, sigma)
            #rho = 0.0004
            init.append(rho)
        for i in range(len(X)):
            init.append(0*init[i])
        for i in range(len(X)):
                #Double the temperature in the middle region
                #Ti = np.exp(-0.1*(Xreduced[i]-Lreduced/2)**2) +1
                if X[i] >= ((L*0.5)-(L*0.02)) and X[i] <= ((L*0.5)+(L*0.02)):
                    Ti = calc_real_T(150,eps)
                    temp.append(Ti)
                else: 
                    Ti = 124
                    temp.append(Ti)
                uThermo = ljs.internal_energy_tv(Ti,hVolume,[moles])
                uThermo = uThermo[0]/moles #J/mol
                uThermo = uThermo/molarMass #J/kg
                #since v=0 we set E = rho*U + 0.5rho*
                init.append(uThermo*rho+0.5*rho*v0**2)
        plt.plot(X, init[2*N:])
        plt.show()

        
        sol1 = odeint(odefunc2FORCELAX, init, tspan)

        solutions = energyMatrixMaker(sol1)
        densities1 = solutions[0]
        velocities1 = solutions[1]
        internalEnergies = solutions[2]
        #pressures1 = pressureIsotermCalc(densities1, T)
        plottingFunctionHeat(densities1, densities1,"rho [kg/m3]", "heatDensity", True, True)
        plottingFunctionHeat(velocities1, densities1,"v [m/s]", "heatVelocity", True, True)
        plottingFunctionHeat(internalEnergies, densities1,"U [J/kg]", "heatInternalEnergy", True, True)




        print("TjoHei!")


    return None


main()

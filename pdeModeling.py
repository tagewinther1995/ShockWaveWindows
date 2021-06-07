import numpy as np
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import os
import imageio
import scipy




def edgeFrontFunction(gridLength, times):
    edgeList = [0]
    for i in range(1,times):
        a = gridLength*i
        edgeList.append(a)
    return edgeList

def edgeBackFunction(gridLength, times):
    edgeList = [N-1]
    for i in range(2,times+1):
        a = gridLength*i
        edgeList.append(a-1)
    return edgeList
def arrayToList(matrix):
    liste = []
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            liste.append(matrix[i][j])
    return liste


N = 200  # number of points to discretize
L = 1.2
xc = 0.3
v = 1
X = np.linspace(0, L, N)  # position along the rod
h = L / (N - 1)
tStop = 0.6
K = 25
tspan = np.linspace(0.0, tStop, K)
k = tStop/(K-1)
realSol = []
dxdt = h/k
dtdx = k/h
for m in range(len(tspan)):
    realSol.append([])
    for n in range(len(X)):
        u = np.exp(-200*(X[n]-xc-v*tspan[m])**2)
        realSol[m].append(u)

npMatrix = np.array(realSol)
solutionList = arrayToList(realSol)


def odefunc(u, t):
    dudt = np.zeros(len(X))

    
    #dudt[-1] = (u[-1] - u[-2])/h|
    #Dss020 differencing scheme
    dudx = []
    #Forward
    #dudt[1] = -v*(-2*u[0]-3*u[1]+6*u[2]-u[3])/(6*h)
    #backward boundary
    #dudt[-2] = -v*(2*u[-1]+3*u[-2]+6*u[-3]+u[-4])/(6*h)
    #Center but removing non existing points

    #dudt[0] = - v*(-u[2] + 8*u[1 ]) / (6*h)
    dudt[1] = - v*(u[2] - u[0]) / (2*h)
    dudt[-2] = - - v*(u[-1] - u[-3]) / (2*h)
    #dudt[-1] = - v*(- 8*u[-2 ] + u[-3]) / (6*h)

    
    


    for i in range(2, N-2):
        #dudt[i] = -v*dudx[i]

        #dudt[i] = - v*(u[i + 1] - u[i - 1]) / (2*h)
        #Testing for fourth order scheme
        dudt[i] = - v*(-u[i + 2] + 8*u[i + 1 ]- 8*u[i - 1 ] + u[i - 2]) / (12*h)
        #dudt[i] = - v*(u[i] - u[i - 1]) / (h)
    dudt[0] = 0
    #dudt[0] =  v*(-u[2] +4* u[1]-3*u[0]) / (2*h) # constant at boundary condition
    #dudt[0] = (u[1] - u[0])/h
    dudt[-1] = 0
    #dudt[-1] = v*(3*u[-1] -4* u[-2]+ u[-3]) / (2*h)
    return dudt


def odefuncLaxFORCE(y, t):
    dudt = np.zeros(len(X))
    dudx = []
    for i in range(1, N-1):
       
        Fplus = 0.5*v*(y[i]+y[i+1])+0.5*dxdt*(y[i]-y[i+1])

        Fminus = 0.5*v*(y[i-1]+y[i])-0.5*dxdt*(y[i]-y[i-1])


        dudt[i] = - (Fplus-Fminus)/h
    dudt[0] = 0
    dudt[-1] = 0
    return dudt

def odefuncFORCE(y, t):
    dudt = np.zeros(len(X))
    dudx = []
    for i in range(1, N-1):
       
        FlfPlus = 0.5*v*(y[i]+y[i+1])+0.5*dxdt*(y[i]-y[i+1])

        FriPlus = v*(0.5*(y[i+1]+y[i])+v*0.5*(1/dxdt)*(y[i]-y[i+1]))
        Fplus = 0.5*(FlfPlus + FriPlus)

        Flf = 0.5*v*(y[i-1]+y[i])-0.5*dxdt*(y[i]-y[i-1])

        Fri = v*(0.5*(y[i-1]+y[i])-v*0.5*(1/dxdt)*(y[i]-y[i-1]))

        Fminus = 0.5*(Flf + Fri)

        dudt[i] = - (Fplus-Fminus)/h
    dudt[0] = 0
    dudt[-1] = 0
    return dudt





def fSolveAdvectionEquation(x, t, u0):
    gridpoints = (len(x) * len(t))
    
    frontList = edgeFrontFunction(N,K)
    backList =  edgeBackFunction(N,K)

    def foo(y):
        tStep = len(u0)
        f = []
        

        for i in range(tStep):
            #initialization setting u0 define u[0] ->u[99]
            fi0 = y[i] - u0[i]
            f.append(fi0)

        for i in range(tStep, gridpoints- tStep):


            #fi = (y[i]-y[i-tStep])/k + v*(y[i]-y[i-1])/h
            
            #y[i-1] og y[i+1] kan være problemer om de overlapper andre grid
            
                

            if i in frontList:
                #Boundary conditions x =0
                
                fi = y[i] - solutionList[i]
            elif i in backList:
                #Boundary conditions x = L
                fi  = y[i] - solutionList[i]
            else:
                #korr
                Fplus = 0.5*v*(y[i]+y[i+1])+0.5*dxdt*(y[i]-y[i+1])
                #korr
                Flf = 0.5*v*(y[i-1]+y[i])-0.5*dxdt*(y[i]-y[i-1])
                #Skrevet inn galt
                #Fri = v*(0.5*(y[i-1]+y[i])-v*0.5*(dxdt)*(y[i]-y[i-1]))

                #skal skrives
                Fri = v*(0.5*(y[i-1]+y[i])-v*0.5*(1/dxdt)*(y[i]-y[i-1]))

                Fminus = 0.5*(Flf + Fri)

                fi = - (y[N+i]+ y[i])- dxdt*(Fplus-Fminus)

            f.append(fi)
        for i in range(gridpoints- tStep,gridpoints ):
            fi = y[i] - solutionList[i]
            f.append(fi)

        #fN = y[len(x)-1] - yN
        # f.append(fN)

        return f
    guess = np.ones(gridpoints)
    yList = scipy.optimize.fsolve(foo, guess)
    #yList = yResult['x'][:]
    yArray = []
    for i in range(len(t)):
        liste = yList[(i*len(t)):(i*len(t)+len(x))]
        yArray.append(liste)

    return yArray


# init = 150.0 * np.ones(X.shape) # initial temperature
# init[0] = 100.0  # one boundary condition
# init[-1] = 200.0 # the other boundary condition
init = []
for i in range(len(X)):
    u0 = np.exp(-200*(X[i]-xc)**2) 
    init.append(u0)

discSol = odeint(odefuncFORCE, init, tspan)

sol = odeint(odefunc, init, tspan)

lax = odeint(odefuncLaxFORCE,init,tspan)


filenames = []

for i in range(0, len(tspan)):
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False
    plt.axis([0, X[-1], (min(realSol[0])+min(realSol[0])*0.1), (max(realSol[0])+max(realSol[0])*0.5)])

    plt.plot(X, sol[i], 'c.', label="MOL: 4.th order Central Differencing")
    plt.plot(X,lax[i],color = 'orange',linestyle="dashed", label = "Lax-Friedrich")
    plt.plot(X,discSol[i], color ="darkgreen",linestyle="dashdot", label = "FORCE")
    plt.plot(X, realSol[i],  color="darkred",linestyle="solid", label= "Analytical")
    plt.xlabel('x [m/s]')
    plt.ylabel('c [mol/m3]')
    plt.legend(loc=2)
    # plt.show()

    stringName = "Python/figures/modelling/" + str(tspan[i]) + 'timeStep.png'
    filenames.append(stringName)
    plt.savefig(stringName)
    plt.clf()

gifName = "Python/figures/modelling/AdvectionMovie.gif"
with imageio.get_writer(gifName, mode='I', duration=0.25) as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)

for m in filenames:
    os.remove(m)

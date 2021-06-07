# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 09:20:16 2021

@author: Tagew
"""
from scipy import constants
from pyctp import ljs_bh
import pandas as pd
import math
import numpy as np
import statistics
import matplotlib.pyplot as plt
from scipy.integrate import simps
import scipy.optimize
from sklearn.linear_model import LinearRegression
import os
import imageio
# Modify system path
import sys
sys.path.append('../pycThermopack/')
# Importing pyThermopack


# Avogadros number
NA = 6.02214076e23


def calc_reduced_T(Ta, eps):
    """ Calculate reduced temperature
    """
    Tstar = np.zeros_like(Ta)
    Tstar = Ta/eps
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

def calc_real_sos(sosStar, m, eps):
    # Calculate reduced speed of sound
    velocity = sosStar/np.sqrt(m/eps)
    return velocity


def calc_Real_P(pStar, sigma, eps):

    pR = np.zeros_like(pStar)
    pR = pStar * eps/(sigma**3)

    return pR


def calc_real_h(hStar, m, eps):
    hR = np.zeros_like(hStar)
    hR = hStar * eps/(m)
    return hR


def calc_reduced_h(hr, m, eps):
    #hStar = np.zeros_like(hr[0])
    h = hr[0]/0.039948  # Change from mole to mass basis
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
    pStar = (pR[0]/eps)*sigma**3
    return pStar


# Her lages en klasse for hver av verdiene der den har et navn en snitt matrise og en standardavvik matrise
class Unit:
    def __init__(self, name, average, deviation, standardError, median):
        self.name = name
        self.average = average
        self.deviation = deviation
        self.standardError = standardError
        self.median = median

# WINDOWSUTGAVE: FOR AT KODEN SKAL FUNKE MÅ ADRESSENE gis med ''\\'', IKKE ''\''


ljs = ljs_bh.ljs_bh()
ljs.init()
ljs.set_tmin(temp=2.0)

# Get parameters for Argon
sigma, eps = ljs.get_sigma_eps()
# For konversjon
aeps = eps*constants.k
m = 6.6335209e-26  # kg


m1 = "/home/tagewm/python/LinuxAnemdFikset/CPROPS1.DAT"
m2 = "/home/tagewm/python/LinuxAnemdFikset/CPROPS2.DAT"
m3 = "/home/tagewm/python/LinuxAnemdFikset/CPROPS3.DAT"
m4 = "/home/tagewm/python/LinuxAnemdFikset/CPROPS4.DAT"
m5 = "/home/tagewm/python/LinuxAnemdFikset/CPROPS5.DAT"
m6 = "/home/tagewm/python/LinuxAnemdFikset/CPROPS6.DAT"
m7 = "/home/tagewm/python/LinuxAnemdFikset/CPROPS7.DAT"
m8 = "/home/tagewm/python/LinuxAnemdFikset/CPROPS8.DAT"
m9 = "/home/tagewm/python/LinuxAnemdFikset/CPROPS9.DAT"
m10 = "/home/tagewm/python/LinuxAnemdFikset/CPROPS10.DAT"

m1 = "/home/tagewm/NEMD/1/CPROPS.DAT"
m2 = "/home/tagewm/NEMD/2/CPROPS.DAT"
m3 = "/home/tagewm/NEMD/3/CPROPS.DAT"
m4 = "/home/tagewm/NEMD/4/CPROPS.DAT"
m5 = "/home/tagewm/NEMD/5/CPROPS.DAT"
m6 = "/home/tagewm/NEMD/6/CPROPS.DAT"
m7 = "/home/tagewm/NEMD/7/CPROPS.DAT"
m8 = "/home/tagewm/NEMD/8/CPROPS.DAT"
m9 = "/home/tagewm/NEMD/9/CPROPS.DAT"
m10 = "/home/tagewm/NEMD/10/CPROPS.DAT"


n1 = "/home/tagewm/python/eqThermoCheck/CPROPS1.DAT"
n2 = "/home/tagewm/python/eqThermoCheck/CPROPS2.DAT"
n3 = "/home/tagewm/python/eqThermoCheck/CPROPS3.DAT"
n4 = "/home/tagewm/python/eqThermoCheck/CPROPS4.DAT"

# Filvariablene fylles inn i filer


files2 = [m1, m2, m3, m4, m5, m6, m7, m8, m9, m10]
eqFiles = [n1, n2, n3, n4]
# Farge rekkefølgen på plottene

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r',
          'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k']

# Navnerekkefølge på dataene du legger inn. (kommer nok til å oppdatere denne delen så man slipper å skrive inn alle variabelene)

names = ['T', 'n', "p", 'x1', 'up', 'h', 'Jm1', 'Jm2',
         'ju1', 'ju2', 'Tx', 'Ty', 'Tz', 'px', 'py', 'pz']
names1 = ['T', 'ρ', "p", 'x1', 'up', 'h', 'Jm1', 'Jm2', 'ju1',
          'ju2', 'Tx', 'Ty', 'Tz', 'px', 'py', 'pz', 'vcm', 'jqflow', 'jq']


# MINOR FUNCTIONS NEEDED TO RUN THE LARGER FUNCTIONS
def timeStamp(array, cLength, time):
    Y = []
    for i in range(len(cLength)):
        Y.append(array[i][time])
    return Y


def regToList(koeff, inter, xList):
    y = []
    for i in range(len(xList)):
        yi = koeff*xList[i] + inter
        y.append(yi)
    return y


def takeClosest(middleValue, y):
    maxValue = max(y)
    maxPos = y.index(maxValue)
    k = 1000
    for i in range(maxPos, len(y)):
        if abs(y[i] - middleValue) < abs(k - middleValue):
            k = y[i]
            index = i
    return index


def findClosestIndex(value, y):
    k = 100000
    for i in range(len(y)):
        if abs(y[i] - value) < abs(k - value):
            k = y[i]
            index = i
    return index


def isNaN(num):
    # sjekker om verdien som leses av er et tall
    return num != num


def columnLength(array):
    # regner ut antall kolonner/lag
    for i in range(len(array)):
        if isNaN(array[i][0]):
            break
    return i


def rowLength(array):
    # lengden på radene, dvs antal tidssteg
    a = len(array[1])-2
    return a


def layerList(array):
    layerList = []
    # henter ut en liste med laglengder
    for i in range(columnLength(array)):
        k = float(array[i][1])
        layerList.append(k)
    return layerList


def noGo(number, variables):
    liste = []
    for i in range(1, variables+1):
        liste.append((i*number)+1)
    return liste


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def intervals(times, length):
    liste = [0, length-1]
    for i in range(1, times):
        liste.append(liste[-1]+1)
        liste.append(liste[-1]+length-1)
    return liste


def additionArray(liste, verdi):
    nyListe = []
    for i in range(len(liste)):
        nyListe.append(liste[i]+verdi)
    return nyListe


def additionMatrix(liste, verdi):
    nyListe = []
    for i in range(len(liste)):
        nyListe.append([])
        for j in range(len(liste[0])):
            nyListe[i].append(liste[i][j]+verdi)
    return nyListe
# END OF MINOR FUNCTIONSLIST

# Funskjonen averageData finner antall lag og regner ut snitt og standardavvik for filene


def averageData(filenames):
    # lager tomme
    values = []
    averages = []
    deviations = []
    errors = []
    medians = []
    # Her hentes ut data fra første filen for å lage en liste over lagene
    data = pd.read_fwf(filenames[0], infer_nrows=1000)
    array0 = data.to_numpy()
    layers = layerList(array0)
    for i in filenames:
        # Her leses datanene fra filene av. De blir gjort om til store matriser og legges inn i en liste values
        data = pd.read_fwf(i, infer_nrows=1000)
        array = data.to_numpy()
        values.append(array)
    for l in range(len(array0)):
        # Her går programmet gjennom fil-matrisene og regner ut snittverdi og standardavvik
        aveRow = []
        sdRow = []
        seRow = []
        medRow = []
        for m in range(2, len(array[l])):
            sample = []
            for k in range(len(filenames)):
                try:
                    if isNaN(float(values[k][l][m])):
                        # sjekker om det
                        break
                    else:
                        b = float(values[k][l][m])
                        sample.append(b)
                except ValueError:
                    print(values[k][l][m])
                    break
            if len(sample) > 0:
                # regner ut snitt og standard avvik
                ave = sum(sample)/len(sample)
                sd = 3*statistics.stdev(sample)
                sa = sd/math.sqrt(len(sample))
                med = statistics.median(sample)
                aveRow.append(ave)
                sdRow.append(sd)
                seRow.append(sa)
                medRow.append(med)
        if len(sdRow) > 0:
            # etter at snitt og SD er regnet ut for raden legges raden inn i snitt og Sd matrise
            averages.append(aveRow)
            deviations.append(sdRow)
            errors.append(seRow)
            medians.append(medRow)
    # Programmet returnerer en liste med lagavstand, en snittmatrise og en sdmatrise
    del aveRow
    del sdRow
    del seRow
    del medRow
    del array
    del data
    return [layers, averages, deviations, errors, medians]


def processing(filenames, names):
    # Formålet med denne funskjonen er å dele snitt og SD matrisen til en liste klasser
    # der hver klasse er en enhet med sitt respektive navn, snitt og sd matrise
    fileRead = averageData(filenames)
    layerList = fileRead[0]
    averages = fileRead[1]
    deviations = fileRead[2]
    errors = fileRead[3]
    medians = fileRead[4]
    values = []
    cLength = len(layerList)
    interval = intervals(len(names), cLength)
    for i in range(len(names)):
        # Her splittes matrisene opp og settes inn i en klasse
        ave = averages[(interval[2*i]):(interval[(2*i)+1]+1)]
        sd = deviations[(interval[2*i]):(interval[(2*i)+1]+1)]
        se = errors[(interval[2*i]):(interval[(2*i)+1]+1)]
        med = medians[(interval[2*i]):(interval[(2*i)+1]+1)]
        name = names[i]
        p = Unit(name, ave, sd, se, med)
        values.append(p)
    # programmet returnerer en liste med lagene og en liste med klassene
    return [layerList, values]


def plottingMD(filenames, names, time):
    # funksjonen tar inn en liste med filer, rekkefølgen på navn og tidspunktet man ønsker å undersøke
    # her hentes dataene ut og plottes
    process = processing(filenames, names)
    layers = process[0]
    values = process[1]
    results = []
    for i in range(len(values)):
        Y = []
        E = []
        ave = values[i].average
        sd = values[i].standardError
        for j in range(len(layers)):
            if len(ave) == len(layers):
                Y.append(ave[j][time])
                E.append(sd[j][time])
            else:
                break
        results.append(Y)
        if len(Y) == len(layers):
            if len(E) == len(layers):
                # Her plottes figurene
                plt.errorbar(layers, Y, E, linestyle='None',
                             marker='.', color=colors[i])
                if min(Y) < 0:
                    plt.axis(
                        [0, layers[-1]+5, (min(Y)+min(Y)*0.1), (max(Y)+max(Y)*0.1)])

                elif max(Y) > 0:
                    plt.axis([0, layers[-1]+5, 0, (max(Y)+max(Y)*0.1)])
                else:
                    plt.axis([0, layers[-1]+5, 0, (max(Y)+0.1)])
                plt.rcParams['axes.spines.right'] = False
                plt.rcParams['axes.spines.top'] = False

                plt.xlabel('x*')
                stringName = "Python/figures/" + values[i].name + '10sim.png'
                yName = values[i].name + '*'
                plt.ylabel(yName)
                # plt.show()
                plt.savefig(stringName, dpi=300)
                plt.clf()
    return None


def writeFileVersion2(filenames, names):
    process = processing(filenames, names)
    layers = process[0]
    values = process[1]
    # legge til lag og liste i alle
    for j in range(len(values)):
        averages = values[j].average
        deviations = values[j].deviation
        standardErrors = values[j].standardError
        medians = values[j].median
        k = 0
        for i in range(len(values[0].average)):
            if k == 64:
                k = 0
            layer = layers[k]
            N = k + 1
            averages[i].insert(0, layer)
            averages[i].insert(0, N)
            deviations[i].insert(0, layer)
            deviations[i].insert(0, N)
            standardErrors[i].insert(0, layer)
            standardErrors[i].insert(0, N)
            medians[i].insert(0, layer)
            medians[i].insert(0, N)
            k = k+1
        values[j].average = averages
        values[j].deviation = deviations
        values[j].standardError = standardErrors
        values[j].median = medians
    aveNames = []
    sdNames = []
    seNames = []
    medNames = []
    for m in range(len(values)):
        aName = values[m].name + 'averages.dat'
        aveNames.append(aName)
        sdName = values[m].name + 'standardDeviatons.dat'
        sdNames.append(sdName)
        seName = values[m].name + 'standardErrors.dat'
        seNames.append(seName)
        medName = values[m].name + 'medians.dat'
        medNames.append(medName)
        head = 'layer \t x* \t' + values[m].name
        np.savetxt(aName, values[m].average, fmt='%.4E', header=head,
                   comments='')
        np.savetxt(sdName, values[m].deviation, fmt='%.4E', header=head,
                   comments='')
        np.savetxt(seName, values[m].standardError, fmt='%.4E', header=head,
                   comments='')
        np.savetxt(medName, values[m].median, fmt='%.4E', header=head,
                   comments='')
    fout = open("Averages.dat", "w")

    for line in open(aveNames[0]):
        fout.write(line)
    for num in range(1, len(aveNames)):
        f = open(aveNames[num])
        for line in f:
            fout.write(line)
        f.close()
    fout.close()
    fout = open("StandardDeviations.dat", "w")
    for line in open(sdNames[0]):
        fout.write(line)
    for num in range(1, len(aveNames)):
        f = open(sdNames[num])
        for line in f:
            fout.write(line)
        f.close()
    fout.close()
    fout = open("StandardErrors.dat", "w")
    for line in open(seNames[0]):
        fout.write(line)
    for num in range(1, len(aveNames)):
        f = open(seNames[num])
        for line in f:
            fout.write(line)
        f.close()
    fout.close()
    fout = open("Medians.dat", "w")
    for line in open(medNames[0]):
        fout.write(line)
    for num in range(1, len(aveNames)):
        f = open(medNames[num])
        for line in f:
            fout.write(line)
        f.close()
    fout.close()

    for m in range(len(values)):
        os.remove(aveNames[m])
        os.remove(sdNames[m])
        os.remove(seNames[m])
        os.remove(medNames[m])
    return None


def findGibbsSurSimp(densities, layers, animation=False, plotting=False):
    # for løkke for relevante tidspunkt
    stop = len(densities[0])-10
    shockPosition = []
    filenames = []
    time = []
    aList = []
    bList = []
    maxYaxis = max(timeStamp(densities, layers, 3))
    for i in range(2, stop):
        # Registrerer start og sluttpunkt
        time.append(i)
        # Henter ut liste for relevant tidspunkt
        data = timeStamp(densities, layers, i)
        # Henter ut makspunktet for lista
        maxValue = max(data)
        # henter ut indexen for makspunktet
        maxPos = data.index(maxValue)
        # Setter et likevektspunkt som er 20 punkter forran makspunktet
        eqPos = maxPos + 15
        # Regner ut snittverdi av bulkfasen forran
        yfAve = (sum(data[-10:-1])/len(data[-10:-1]))
        yAveList = []
        # lager en snittverdi liste for bulkfasen foran
        for j in range(len(layers)):
            yAveList.append(yfAve)
        # Finner middelverdipunktet
        cV = (maxValue + data[eqPos])/2
        # Finner indexen nærmest middelverdipunktet
        c = takeClosest(cV, data)
        # Setter et punkt a og b forann og bak C punktet
        a = c - 5
        aList.append(a)
        b = c + 5
        bList.append(b)

        aN = 10
        bN = 10
        # Gjør en simpsonintegrasjon
        I1 = simps(data[a-aN:b+bN+1], layers[a-aN:b+bN+1])
        # Making linear regression lines for the front(f) and back(b)
        # kanskje lurere å basere frontlinjen på et snitt av likevektsverdiene
        yb = np.array(data[(a-aN):a+1])
        xb = np.array(layers[(a-aN):a+1]).reshape((-1, 1))
        #yf = np.array(data[b:(b+bN+1)])
        #xf = np.array(layers[b:(b+bN+1)]).reshape((-1, 1))

        # Create a model and fit it
        modelBack = LinearRegression()
        modelBack.fit(xb, yb)
        #modelFront = LinearRegression()
        # modelFront.fit(xf,yf)

        aBack = modelBack.coef_
        bBack = modelBack.intercept_
        #aFront = modelFront.coef_
        #bFront = modelFront.intercept_
        yBack = regToList(aBack, bBack, layers)
        #yFront= regToList(aFront,bFront,layers)

        def f(l):

            f = (I1-(0.5*aBack*(l*l-layers[a-aN]**2) + bBack *
                     (l - layers[a-aN]) + yfAve*(layers[b+bN] - l)))
            #d = l - layers[a-aN]
            #f = (I1 -(0.5*float(aBack)*(l*l-(l-d)**2) + float(bBack)*(l - (l-d)) + yfAve*((l+d) -l)))
            return f

        def m(l):
            # Bjørn:(((l-layers[a-aN])*(yBack[a-aN]+yBack[a])/2+(layers[(b+bN)]-l)*(yfAve)-I1)**2)
            #((0.5*aBack*(l*l-layers[a-aN]**2) + bBack*(l - layers[a-aN]) + yfAve*(layers[b+bN] -l) - I1)**2)
            return (((0.5*aBack*(l*l-layers[a-aN]**2) + bBack*(l - layers[a-aN]) + yfAve*(layers[b+bN] - l) - I1))**2)

        guess = layers[maxPos]
        # Root method
        solution = scipy.optimize.root(f, guess)
        # Minimize method
        #solMin = scipy.optimize.minimize_scalar(m, bounds=(0, layers[-1]),method='bounded')

        shockPosition.append(solution['x'][0])
        if plotting == True:
            #line = solution['x'][0]
            line2 = solution['x'][0]
            plt.rcParams['axes.spines.right'] = False
            plt.rcParams['axes.spines.top'] = False
            plt.xlabel('x*')
            plt.ylabel('ρ*')
            plt.plot(layers, yBack, 'b')
            plt.plot(layers, yAveList, 'r')
            plt.plot(layers, data, 'g.', label='Density')
            plt.axvline(x=line2)
            plt.axvline(x=layers[a-10])
            plt.axvline(x=layers[b+10])

            plt.axis([0, layers[-1], 0, (maxYaxis+maxYaxis*0.25)])
            stringName = 'Python/figures/gibbsEqSur/' + str(i) + 'timeStep.png'
            filenames.append(stringName)
            plt.savefig(stringName, dpi=300)
            plt.clf()
            #l = solution['x'][0]

    if animation == True and plotting == True:
        with imageio.get_writer('Python/figures/gibbsEqSur/movie.gif', mode='I', duration=1) as writer:
            for filename in filenames:
                image = imageio.imread(filename)
                writer.append_data(image)

    return [shockPosition, time, aList, bList]


def shockSpeed(timeLayers, shockPosition, degree, plot=False):
    # This function returns the shock velocity for a given set of shockwave positions and timesteps
    poly = np.polyfit(timeLayers, shockPosition, degree)
    ployDer = np.polyder(poly)
    shockVelocity = []
    polyShock = []
    for i in range(len(timeLayers)):
        a = 0
        b = 0
        for j in range(len(ployDer)):
            a = a + ployDer[j]*timeLayers[i]**(degree-j-1)

        shockVelocity.append(a)
        if plot == True:
            for m in range(len(poly)):
                b = b + poly[m]*timeLayers[i]**(degree-m)
            polyShock.append(b)
    if plot == True:
        plt.rcParams['axes.spines.right'] = False
        plt.rcParams['axes.spines.top'] = False
        plt.xlabel('t*')
        plt.ylabel('l*')
        plt.plot(timeLayers, shockPosition, 'b.')
        plt.plot(timeLayers, polyShock, 'r')
        plt.savefig("shockwaveSpeedfit.PNG", dpi=300)
        plt.show()
        plt.clf()
    return shockVelocity


def linearRegression(xList, yList, filename, xAxisName, yAxisName, linReg=True, farge="c."):
    yb = np.array(yList)
    xb = np.array(xList).reshape((-1, 1))
    modelBack = LinearRegression()
    modelBack.fit(xb, yb)
    aBack = modelBack.coef_
    bBack = modelBack.intercept_
    yBack = regToList(aBack, bBack, xList)
    plt.ylabel(yAxisName)
    plt.xlabel(xAxisName)
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False
    #plt.axis([0,xList[-1]+ xList[-1]*0.1,0,(max(yList)+max(yList)*0.25)])
    plt.plot(xList, yList, farge)
    if linReg == True:
        plt.plot(xList, yBack)
    plt.savefig(filename, dpi=300)
    plt.show()
    return None

def linearRegressionError(xList, yList, errorList, filename, xAxisName, yAxisName, linReg=True):
    yb = np.array(yList)
    xb = np.array(xList).reshape((-1, 1))
    modelBack = LinearRegression()
    modelBack.fit(xb, yb)
    aBack = modelBack.coef_
    bBack = modelBack.intercept_
    yBack = regToList(aBack, bBack, xList)
    plt.ylabel(yAxisName)
    plt.xlabel(xAxisName)
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False
    #plt.axis([0,xList[-1]+ xList[-1]*0.1,0,(max(yList)+max(yList)*0.25)])
    plt.errorbar(xList, yList, errorList, linestyle='None',marker='.', color= 'c')
    if linReg == True:
        plt.plot(xList, yBack)
    plt.savefig(filename, dpi=300)
    plt.show()
    return None

def animateArray(classType, layers, shockPosition):
    array = classType.average
    se = classType.standardError
    maxYaxis = max(timeStamp(array, layers, 4))

    filenames = []
    for i in range(len(array[0])):
        yList = timeStamp(array, layers, i)
        eList = timeStamp(se, layers, i)
        plt.rcParams['axes.spines.right'] = False
        plt.rcParams['axes.spines.top'] = False
        plt.xlabel('x*')
        plt.ylabel(classType.name + '*')
        plt.errorbar(layers, yList, eList, linestyle='None', marker='.')
        if len(shockPosition) > 0:
            if i >= 4:
                if i <= len(array[0])-4:
                    try:
                        plt.axvline(x=shockPosition[i-4])
                    except IndexError:
                        print("øut of range")

        plt.axis([0, layers[-1]+0.1*layers[-1], 0, (maxYaxis+maxYaxis*0.25)])
        stringName = classType.name + str(i) + 'timeStep.png'
        filenames.append(stringName)
        plt.savefig(stringName, dpi=300)
        plt.clf()
    gifName = "Python/figures/animate/" + classType.name + 'movie.gif'
    with imageio.get_writer(gifName, mode='I', duration=0.5) as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

    for m in filenames:
        os.remove(m)
    return None


def animateArrayEOScomparison(classType, eosArray, layers, shockPosition):
    array = classType.average
    se = classType.standardError
    maxYaxis = max(timeStamp(array, layers, 4))

    filenames = []
    for i in range(len(array[0])):
        yList = timeStamp(array, layers, i)
        #yEOS = timeStamp()
        eList = timeStamp(se, layers, i)
        plt.rcParams['axes.spines.right'] = False
        plt.rcParams['axes.spines.top'] = False
        plt.xlabel('x*')
        plt.ylabel(classType.name + '*')
        plt.errorbar(layers, yList, eList, linestyle='None', marker='.')
        if len(shockPosition) > 0:
            if i >= 4:
                if i <= len(array[0])-4:
                    try:
                        plt.axvline(x=shockPosition[i-4])
                    except IndexError:
                        print("øut of range")

        plt.axis([0, layers[-1]+0.1*layers[-1], 0, (maxYaxis+maxYaxis*0.25)])
        stringName = classType.name + str(i) + 'timeStep.png'
        filenames.append(stringName)
        plt.savefig(stringName, dpi=300)
        plt.clf()
    gifName = "Python/figures/animate/" + classType.name + 'movie.gif'
    with imageio.get_writer(gifName, mode='I', duration=0.5) as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

    for m in filenames:
        os.remove(m)
    return None


def speedOfSoundPressure(T, pressureList, phase):
    ljs = ljs_bh.ljs_bh()
    ljs.init()
    ljs.set_tmin(temp=2.0)
    z = [1]
    # Get parameters for Argon
    sigma, eps = ljs.get_sigma_eps()
    # For konversjon
    aeps = eps*constants.k
    m = 6.6335209e-26
    sosList = []
    for i in range(len(pressureList)):
        TeqR = calc_real_T(T, eps)
        Peq = pressureList[i]
        PeqR = calc_Real_P(Peq, sigma, aeps)
        if phase == ljs.LIQPH:
            sos = ljs.speed_of_sound(TeqR, PeqR, z, z, z, 0, 1, ljs.LIQPH)
        else:
            sos = ljs.speed_of_sound(TeqR, PeqR, z, z, z, 1, 0, ljs.VAPPH)
        sosStar = calc_reduced_sos(sos, m, aeps)
        sosList.append(sosStar)
    return sosList


def linearRegressionComparison(xList, yListEOS, yListNEMD, filename, xAxisName, yAxisName):
    # in case of errorbars
    #ex = [10*i for i in ex]
    #ey = [10*i for i in ex]

    yb1 = np.array(yListEOS)
    yb2 = np.array(yListNEMD)
    xb = np.array(xList).reshape((-1, 1))
    modelBack = LinearRegression()
    modelBack.fit(xb, yb1)
    aBack1 = modelBack.coef_
    bBack1 = modelBack.intercept_
    yBack1 = regToList(aBack1, bBack1, xList)
    modelBack = LinearRegression()
    modelBack.fit(xb, yb2)
    aBack2 = modelBack.coef_
    bBack2 = modelBack.intercept_
    yBack2 = regToList(aBack2, bBack2, xList)
    plt.ylabel(yAxisName)
    plt.xlabel(xAxisName)
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False
    #plt.axis([0,xList[-1]+ xList[-1]*0.1,0,(max(yList)+max(yList)*0.25)])
    #plt.errorbar(xList, yList1, yerr = ey,xerr =ex , linestyle='None', marker='.',color="c",label= "NEMD simulations")

    plt.plot(xList, yListEOS, 'gx', label="3 Order pertubation")
    plt.plot(xList, yListNEMD, 'c.', label="Nemd data")
    

    #plt.plot(xList, yBack1,"c")
    #plt.plot(xList,yBack2, "g")
    plt.legend()
    plt.savefig(filename, dpi=300)
    plt.show()
    return None
def linearTempComparison(xList, yListTx, yListTy,yListTz, filename, xAxisName, yAxisName):
    # in case of errorbars
    #ex = [10*i for i in ex]
    #ey = [10*i for i in ex]

    plt.ylabel(yAxisName)
    plt.xlabel(xAxisName)
    plt.axis([0, layers[-1]+5, 0, (max(yListTx)+max(yListTx)*0.1)])
    #plt.axis([0,xList[-1]+ xList[-1]*0.1,0,(max(yList)+max(yList)*0.25)])
    #plt.errorbar(xList, yList1, yerr = ey,xerr =ex , linestyle='None', marker='.',color="c",label= "NEMD simulations")

    
    plt.plot(xList, yListTy, 'cx', label="Tyy")
    plt.plot(xList, yListTz, 'r+', label="Tzz")
    plt.plot(xList, yListTx, 'g.', label="Txx")
    
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False
    #plt.plot(xList, yBack1,"c")
    #plt.plot(xList,yBack2, "g")
    plt.legend()
    plt.savefig(filename, dpi=300)
    plt.show()
    return None

def linearRegressionComparisonError(xList, yListEOS, yListNEMD, errorList, filename, xAxisName, yAxisName):
    # in case of errorbars
    #ex = [10*i for i in ex]
    #ey = [10*i for i in ex]

    yb1 = np.array(yListEOS)
    yb2 = np.array(yListNEMD)
    xb = np.array(xList).reshape((-1, 1))
    modelBack = LinearRegression()
    modelBack.fit(xb, yb1)
    aBack1 = modelBack.coef_
    bBack1 = modelBack.intercept_
    yBack1 = regToList(aBack1, bBack1, xList)
    modelBack = LinearRegression()
    modelBack.fit(xb, yb2)
    aBack2 = modelBack.coef_
    bBack2 = modelBack.intercept_
    yBack2 = regToList(aBack2, bBack2, xList)
    plt.ylabel(yAxisName)
    plt.xlabel(xAxisName)
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False
    #plt.axis([0,xList[-1]+ xList[-1]*0.1,0,(max(yList)+max(yList)*0.25)])
    #plt.errorbar(xList, yList1, yerr = ey,xerr =ex , linestyle='None', marker='.',color="c",label= "NEMD simulations")

    plt.plot(xList, yListEOS, 'gx', label="3 Order pertubation")
    #plt.plot(xList, yListNEMD, 'c.', label="Nemd data")
    plt.errorbar(xList, yListNEMD, errorList, linestyle='None',marker='.', color="c",label="Nemd data")

    #plt.plot(xList, yBack1,"c")
    #plt.plot(xList,yBack2, "g")
    plt.legend()
    plt.savefig(filename, dpi=300)
    plt.show()
    return None


def excessDensity(layers, densityProperty, densities, shockPosition, intStart, intStopp, animation=False, plotting=False):
    # for løkke for relevante tidspunkt
    stop = len(densityProperty.average[0])-10
    excessProperty = []
    filenames = []
    time = []
    excessDens = []
    maxYaxis = max(timeStamp(densityProperty.average, layers, 5))
    minYaxis = min(timeStamp(densityProperty.average, layers, 5))
    for k in range(len(densities)):
        excessDens.append([])
        for l in range(len(densities[0])):
            dens = densityProperty.average[k][l] * densities[k][l]
            excessDens[k].append(dens)

    # Legger til variabel for å ha kontroll på sjokk posisjon
    k = 0
    for i in range(2, stop):
        time.append(i)
        l = shockPosition[k]
        data = timeStamp(excessDens, layers, i)
        maxValue = max(data)
        # bestemmer likevekts system utfra sjokkposisjon
        eqPos = l + 15
        yfAve = (sum(data[-10:-1])/len(data[-10:-1]))
        yAveList = []
        for j in range(len(layers)):
            yAveList.append(yfAve)

        #c = findClosestIndex(l,layers)
        a = intStart[k]
        b = intStopp[k]
        aN = 10
        bN = 10
        I1 = simps(data[a-aN:b+bN+1], layers[a-aN:b+bN+1])
        # Making linear regression lines for the front(f) and back(b)
        # kanskje lurere å basere frontlinjen på et snitt av likevektsverdiene
        yb = np.array(data[(a-aN):a+1])
        xb = np.array(layers[(a-aN):a+1]).reshape((-1, 1))
        #yf = np.array(data[b:(b+bN+1)])
        #xf = np.array(layers[b:(b+bN+1)]).reshape((-1, 1))
        # Create a model and fit it
        modelBack = LinearRegression()
        modelBack.fit(xb, yb)
        #modelFront = LinearRegression()
        # modelFront.fit(xf,yf)
        aBack = modelBack.coef_
        bBack = modelBack.intercept_
        rNumber = modelBack.score(xb, yb)
        # print(rNumber)
        #aFront = modelFront.coef_
        #bFront = modelFront.intercept_
        yBack = regToList(aBack, bBack, layers)
        #yFront= regToList(aFront,bFront,layers)
        # To get symmetric integration, introduce a value of d
        d = l - layers[a-aN]
        excessValue = I1 - (0.5*float(aBack)*(l*l-layers[a-aN]**2) + float(
            bBack)*(l - layers[a-aN]) + yfAve*(layers[b+bN] - l))
        #excessValue = I1 -(0.5*float(aBack)*(l*l-(l-d)**2) + float(bBack)*(l - (l-d)) + yfAve*((l+d) -l))
        excessProperty.append(excessValue)
        k = k + 1
        if plotting == True:
            #line = solution['x'][0]
            line2 = l
            plt.rcParams['axes.spines.right'] = False
            plt.rcParams['axes.spines.top'] = False
            plt.xlabel("x*")
            plt.ylabel("ρ" + densityProperty.name + '*')
            plt.plot(layers, yBack, 'b')
            plt.plot(layers, yAveList, 'r')
            plt.plot(layers, data, 'g.', label='Density')
            plt.axvline(x=line2)
            plt.axvline(x=layers[a-10])
            plt.axvline(x=layers[b+10])
            if minYaxis < 0:
                plt.axis(
                    [0, layers[-1], (minYaxis - abs(minYaxis*0.25)), 0.5*(maxYaxis)])
            else:
                plt.axis([0, layers[-1], 0, (maxYaxis+maxYaxis*0.25)])
            stringName = "Python/figures/excessProperty/" + \
                str(i) + 'timeStep.png'
            filenames.append(stringName)
            plt.savefig(stringName, dpi=300)
            plt.clf()

    if animation == True and plotting == True:
        with imageio.get_writer('Python/figures/excessProperty/movie.gif', mode='I', duration=1) as writer:
            for filename in filenames:
                image = imageio.imread(filename)
                writer.append_data(image)

    return [excessProperty, time]


def EOSentropyDerivation(temperatures, pressures):
    ljs = ljs_bh.ljs_bh()
    ljs.init()
    ljs.set_tmin(temp=2.0)

    # Get parameters for Argon
    sigma, eps = ljs.get_sigma_eps()
    # For konversjon
    aeps = eps*constants.k
    m = 6.6335209e-26  # kg
    # Plot phase envelope using a1-a3
    z = np.array([1])
    entropyArray = []
    for i in range(len(temperatures)):
        entropyArray.append([])
        for j in range(len(temperatures[0])):
            realT = calc_real_T(temperatures[i][j], eps)
            realP = calc_Real_P(pressures[i][j], sigma, aeps)
            entropyReal = ljs.entropy(realT, realP, z, ljs.LIQPH)
            entropyStar = calc_reduced_s(entropyReal, m)
            entropyArray[i].append(entropyStar)
    return entropyArray


def EOSenthalpyDerivation(temperatures, pressures):
    ljs = ljs_bh.ljs_bh()
    ljs.init()
    ljs.set_tmin(temp=2.0)

    # Get parameters for Argon
    sigma, eps = ljs.get_sigma_eps()
    # For konversjon
    aeps = eps*constants.k
    m = 6.6335209e-26  # kg
    z = np.array([1])
    enthalpyArray = []
    for i in range(len(temperatures)):
        enthalpyArray.append([])
        for j in range(len(temperatures[0])):
            realT = calc_real_T(temperatures[i][j], eps)
            realP = calc_Real_P(pressures[i][j], sigma, aeps)
            enthalpyReal = ljs.enthalpy(realT, realP, z, ljs.LIQPH)
            enthalpyStar = calc_reduced_h(enthalpyReal, m, aeps)
            enthalpyArray[i].append(enthalpyStar)

    return enthalpyArray


def EOSinternalEnergyDerivation(temperatures, densities, Np):
    ljs = ljs_bh.ljs_bh()
    ljs.init()
    ljs.set_tmin(temp=2.0)

    # Get parameters for Argon
    sigma, eps = ljs.get_sigma_eps()
    # For konversjon
    aeps = eps*constants.k
    #density = eqDensity
    # Plot phase envelope using a1-a3
    #moles = Np/NA
    m = 6.6335209e-26  # kg
    #volumeStar = Np/density
    volumeStar = 2.353*4.668*4.668
    #volumeReal = calc_Real_volume(volumeStar,sigma)
    internalEnergyArray = []
    for i in range(len(temperatures)):
        internalEnergyArray.append([])
        for j in range(len(temperatures[0])):
            np = densities[i][j]*volumeStar
            moles = np/NA
            #volumeStar = Np/densities[i][j]
            volumeReal = calc_Real_volume(volumeStar, sigma)
            realT = calc_real_T(temperatures[i][j], eps)
            Ur = ljs.internal_energy_tv(realT, volumeReal, [moles])
            uStar = calc_reduced_u(Ur, m, aeps, np)
            internalEnergyArray[i].append(uStar)
    return internalEnergyArray


def EOSpressureDerivation(temperatures, densities, Np):
    ljs = ljs_bh.ljs_bh()
    ljs.init()
    ljs.set_tmin(temp=2.0)

    # Get parameters for Argon
    sigma, eps = ljs.get_sigma_eps()
    # For konversjon
    aeps = eps*constants.k
    #density = eqDensity
    volumeStar = 2.353*4.668*4.668
    #moles = Np/NA
    m = 6.6335209e-26  # kg
    #volumeStar = Np/density
    #volumeReal = calc_Real_volume(volumeStar,sigma)
    pressureArray = []
    for i in range(len(temperatures)):
        pressureArray.append([])
        for j in range(len(temperatures[0])):
            moles = densities[i][j]*volumeStar/NA

            #volumeStar = Np/densities[i][j]
            volumeReal = calc_Real_volume(volumeStar, sigma)
            realT = calc_real_T(temperatures[i][j], eps)
            pReal = ljs.pressure_tv(realT, volumeReal, [moles])
            pStar = calc_reduced_p(pReal, sigma, aeps)
            pressureArray[i].append(pStar)
    return pressureArray


def surfaceTemperature(excessEntropies, excessInternalEnergies, times, degree, plot=False):
    # This function returns the shock velocity for a given set of shockwave positions and timesteps
    poly = np.polyfit(excessEntropies, excessInternalEnergies, degree)
    ployDer = np.polyder(poly)
    surfaceTemp = []
    polyShock = []
    for i in range(len(excessEntropies)):
        a = 0
        b = 0
        for j in range(len(ployDer)):
            a = a + ployDer[j]*excessEntropies[i]**(degree-j-1)

        surfaceTemp.append(a)
        if plot == True:
            for m in range(len(poly)):
                b = b + poly[m]*excessEntropies[i]**(degree-m)
            polyShock.append(b)
    if plot == True:
        plt.rcParams['axes.spines.right'] = False
        plt.rcParams['axes.spines.top'] = False

        plt.xlabel("$ρ_s^s*$")
        plt.ylabel("$ρ_u^s*$")
        plt.plot(excessEntropies, excessInternalEnergies, 'b.')
        plt.plot(excessEntropies, polyShock, 'r')
        plt.savefig("surfaceTempDer.PNG", dpi=300)
        plt.show()
        plt.clf()

    return surfaceTemp


def nemdInternalEnergy(up, temp):
    u = []
    for i in range(len(temp)):
        u.append([])
        for j in range(len(temp[0])):
            uInternal = up[i][j] + 1.5 * temp[i][j]
            u[i].append(uInternal)
    return u


def nemdEOSdifferences(EOSList, nemdList):
    differences = []
    for i in range(len(EOSList)):
        diff = EOSList[i]-nemdList[i]
        differences.append(diff)
    return differences


def specificEntropyCalc(filename):
    data = pd.read_fwf(filename, infer_nrows=1000)
    array = data.to_numpy()
    b = array.tolist()
    Np = 16000
    fix = []

    for i in range(len(array)):
        l = b[i][0].split('\t')
        fix.append(l)
    for i in range(len(fix)):
        for j in range(len(fix[0])):
            fix[i][j] = float(fix[i][j])

    pressures = timeStamp(fix, array, 3)
    densities = timeStamp(fix, array, 2)
    temperatures = timeStamp(fix, array, 1)
    ljs = ljs_bh.ljs_bh()
    ljs.init()
    ljs.set_tmin(temp=2.0)
    sigma, eps = ljs.get_sigma_eps()
    aeps = eps*constants.k
    m = 6.6335209e-26
    z = np.array([1])
    entropyArray = []
    entropyArrayTV = []
    moles = Np/NA

    for i in range(len(temperatures)):
        volumeStar = Np/densities[i]
        volumeReal = calc_Real_volume(volumeStar, sigma)
        realT = calc_real_T(temperatures[i], eps)
        realP = calc_Real_P(pressures[i], sigma, aeps)
        entropyReal = ljs.entropy(realT, realP, z, ljs.LIQPH)
        entropyRealTV = ljs.entropy_tv(realT, volumeReal, [moles])
        print(entropyRealTV)
        entropyStar = calc_reduced_s(entropyReal, m)
        entropyArray.append(entropyStar)
    pressureArray = []
    for i in range(len(temperatures)):
        volumeStar = Np/densities[i]
        volumeReal = calc_Real_volume(volumeStar, sigma)
        realT = calc_real_T(temperatures[i], eps)
        pReal = ljs.pressure_tv(realT, volumeReal, [moles])
        pStar = calc_reduced_p(pReal, sigma, aeps)
        pressureArray.append(pStar)
    linearRegression(temperatures, entropyArray,
                     "specificEntropyBjorn.png", "T*", "s*", False, "c.")
    linearRegressionComparison(
        densities, pressureArray, pressures, "pressureComparisonBjorn.png", "n*", "p*")
    # lagrer entropiene
    for i in range(len(fix)):
        fix[i].append(entropyArray[i])
        fix[i].append(pressureArray[i])
    head = 'case \t T* \t rho* \t p* \t s \t pEOS'
    np.savetxt("termopackEntropies.DAT", fix,
               fmt='%.4E', header=head, comments='')

    return [entropyArray, pressureArray]


def machNumberCalc(timeLayers, shockPosition, eqTemp, eqPress, plot=False):
    time = []
    for i in range(len(timeLayers)):
        time.append(timeLayers[i]+0.5)
    #shockVelocity = shockSpeed(timeLayers, shockPosition, 3)
    shockVelocity = shockSpeed(timeLayers, shockPosition, 3)
    ljs = ljs_bh.ljs_bh()
    ljs.init()
    ljs.set_tmin(temp=2.0)
    sigma, eps = ljs.get_sigma_eps()
    aeps = eps*constants.k
    m = 6.6335209e-26
    TeqR = calc_real_T(eqTemp, eps)
    PeqR = calc_Real_P(eqPress, sigma, aeps)
    z = [1]
    sos = ljs.speed_of_sound(TeqR, PeqR, z, z, z, 0, 1, ljs.LIQPH)
    sosStar = calc_reduced_sos(sos, m, aeps)
    machNumbers = []
    for i in range(len(shockVelocity)):
        ma = shockVelocity[i]/sosStar
        machNumbers.append(ma)
    if plot == True:
        #linearRegression(timeLayers, machNumbers,"machNumbers.png", "t*", "M", False, "k.")
        linearRegression(time, machNumbers,"machNumbers.png", "t*", "M", False, "k.")
    return [machNumbers, sosStar]


def sackurTetrodeEntropy(temperature, density):
    ljs = ljs_bh.ljs_bh()
    ljs.init()
    ljs.set_tmin(temp=2.0)

    # Get parameters for Argon
    sigma, eps = ljs.get_sigma_eps()
    # For konversjon
    aeps = eps*constants.k
    m = 6.6335209e-26  # kg
    # Plot phase envelope using a1-a3
    entropyArray = []
    for i in range(len(temperature)):
        entropyArray.append([])
        for j in range(len(temperature[0])):
            realT = calc_real_T(temperature[i][j], eps)
            #lamb = constants.Planck / math.sqrt(2*math.pi*m*constants.k*realT)
            #entropyReal = constants.k*(5/2 - math.log(density[i][j]*lamb**3/(sigma**3)))
            entropyStar = (
                3*math.log(temperature[i][j])/2 - math.log(density[i][j]))
            entropyArray[i].append(entropyStar)
    return entropyArray

# SETT HVILKET TIDSPUNKT DU VIL KJØRE
tidspunkt = 10
# Skriver snitt, SD, SE og Middelverdi til DAT.filer

#writeFileVersion1(filer, names)
# writeFileVersion2(files2,names1)
#fil = "termopak.txt"
#array = specificEntropyCalc(fil)

# lager et sett med klasser
plottingMD(files2,names1,tidspunkt)
#writeFileVersion2(files2,names1)
process = processing(files2, names1)
layers = process[0]
units = process[1]
density = units[1].average



# Bestemmer sjokkbølge posisjon
shockFront = findGibbsSurSimp(density, layers, True, True)
shockPosition = shockFront[0]
timeList = shockFront[1]
time = []
for i in range(len(timeList)):
    time.append(timeList[i]+0.5)
#b = shockSpeed(shockFront[1], shockFront[0], 3, True)
px = timeStamp(units[13].average,layers,tidspunkt)
py = timeStamp(units[14].average,layers,tidspunkt)
pz = timeStamp(units[15].average,layers,tidspunkt)
p = timeStamp(units[2].average,layers,tidspunkt)
pAV = []
for i in range(len(p)):
    pSnitt = (px[i]+py[i]+pz[i])/3
    pAV.append(pSnitt)

linearRegressionComparison(layers,p,pAV,"pressureTest.png", "x*","p")


b = shockSpeed(time, shockFront[0], 3, True)
intStart = shockFront[2]
intStopp = shockFront[3]

# Beregner lydhastighet og mach tall
eqTemp = 1
pressures = units[2].average
p = timeStamp(pressures, layers, 1)
eqPress = sum(p[-15:])/len(p[-15:])
machSpeedOfSound = machNumberCalc(
    timeList, shockPosition, eqTemp, eqPress, True)

#temperature check
Tx = timeStamp(units[10].average,layers,tidspunkt)
Ty = timeStamp(units[11].average,layers,tidspunkt)
Tz = timeStamp(units[12].average,layers,tidspunkt)
linearTempComparison(layers,Tx,Ty,Tz,"Tempcompar.png", "x*", "T*")

# Beregning av entropi ved bruk av thermopakke EOS og Utregning av excess entropi
temperatures = units[0].average

entropies = EOSentropyDerivation(temperatures, pressures)
entropiesIdeal = sackurTetrodeEntropy(temperatures, density)
sEOS = timeStamp(entropies, layers, 1)
eqS = sum(sEOS[-15:])/len(sEOS[-15:])
sSack = timeStamp(entropiesIdeal, layers, 1)
eqSack = sum(sSack[-15:])/len(sSack[-15:])
diff = eqSack - eqS
ADentropies = additionMatrix(entropies, diff)


f = Unit("s", entropies, [], [], [])
units.append(f)
excessEntropy = excessDensity(
    layers, units[-1], density, shockFront[0], intStart, intStopp,True, True)
excessTime = []
for i in range(len(excessEntropy[1])):
    excessTime.append(excessEntropy[1][i]+0.5)
#linearRegression(excessEntropy[1], excessEntropy[0],"excessEntropyDensity.png", "t*", "$ρ_s^s*$", False, "r.")
linearRegression(excessTime, excessEntropy[0],"excessEntropyDensity.png", "t*", "$ρ_s^s*$", False, "r.")
entropyEOSPoint = timeStamp(entropies, layers, tidspunkt)
#entropyEOSPoint = additionArray(entropyEOSPoint,diff)
entropySackurPoint = timeStamp(entropiesIdeal, layers, tidspunkt)
linearRegressionComparison(
    layers, entropyEOSPoint, entropySackurPoint, "entropyShockComparison.png", "x*", "s*")
linearRegression(layers, entropyEOSPoint, "EOSentropy.png",
                 "x*", "s*", False, "r.")


# Sammenligning av NEMD og EOS data: Entalpi
enthalpies = EOSenthalpyDerivation(temperatures, pressures)
nemdEnthalpy = timeStamp(units[5].average, layers, tidspunkt)
nemdEnthalpyError = timeStamp(units[5].standardError, layers, tidspunkt)
enthalpiPoint = timeStamp(enthalpies, layers, tidspunkt)

linearRegressionComparisonError(layers, enthalpiPoint,
                           nemdEnthalpy,nemdEnthalpyError, "enthalpyComparison.png", "x*", "h*")
comparison = nemdEOSdifferences(enthalpiPoint, nemdEnthalpy)
linearRegression(layers, comparison,
                 "enthalpyShockComparisonDifferences.png", "x*", "Δh*", False, "c.")


# Sammenligning av NEMD og EOS data: Indre energi:

u = EOSinternalEnergyDerivation(temperatures, density, 16000)
nemdInternal = nemdInternalEnergy(units[4].average, temperatures)
nemdInternalError = nemdInternalEnergy(units[4].standardError, units[0].standardError)
nemdInteralPoint = timeStamp(nemdInternal, layers, tidspunkt)
nemdInternalErrorPoint = timeStamp(nemdInternalError, layers, tidspunkt)
uPoint = timeStamp(u, layers, tidspunkt)
linearRegressionComparisonError(layers, uPoint, nemdInteralPoint,nemdInternalErrorPoint,
                           "InternalEnergyShockComparison.png", "x*", "u*")
comparison = nemdEOSdifferences(uPoint, nemdInteralPoint)
linearRegression(layers, comparison,
                 "InternalEnergyShockComparisonDifferences.png", "x*", "Δu*", False, "c.")
f = Unit("u", nemdInternal, [], [], [])
excessU = excessDensity(
    layers, f, density, shockFront[0], intStart, intStopp,True,True)
#linearRegression(excessU[1], excessU[0],"excessInternalEnergy.png", "t*", "$ρ_u^s*$", False, "m.")
linearRegression(excessTime, excessU[0],"excessInternalEnergy.png", "t*", "$ρ_u^s*$", False, "m.")

#ENTALPI DOBBELSJEKK
entalpiSjekk = []

densPoint = timeStamp(units[1].average,layers,tidspunkt)
pPoint = timeStamp(units[2].average,layers,tidspunkt)
for i in range(len(nemdInteralPoint)):
    h = nemdInteralPoint[i] + pPoint[i]/densPoint[i]
    entalpiSjekk.append(h)
linearRegressionComparison(layers,enthalpiPoint,entalpiSjekk,"entalpiTest.png", "x*","h")


# Sammenligning av EOS of nemd trykk
pEOS = EOSpressureDerivation(temperatures, density, 16000)
pEOSpoint = timeStamp(pEOS, layers, tidspunkt)
pNEMD = timeStamp(units[2].average, layers, tidspunkt)
pNEMDerror = timeStamp(units[2].standardError, layers, tidspunkt)
linearRegressionComparison(layers, pEOSpoint, pNEMD,
                           "pressureShockComparison.png", "x*", "p*")

linearRegressionComparisonError(layers, pEOSpoint, pNEMD, pNEMDerror,
                           "pressureShockComparison.png", "x*", "p*")
                           
comparison = nemdEOSdifferences(pEOSpoint, pNEMD)
linearRegression(layers, comparison,
                 "pressureShockComparisonDifferences.png", "x*", "Δp*", False)

# Determining surface Temperature

surfaceT = surfaceTemperature(excessEntropy[0][-5:], excessU[0][-5:], excessEntropy[1][-5:], 1, True)

#average density and pressure for utilizing uv flash
internalArray = timeStamp(nemdInternal, layers, 1)
eqU = sum(internalArray[-15:])/len(internalArray[-15:])
densArray = timeStamp(density, layers, 1)
eqn = sum(densArray[-15:])/len(densArray[-15:])
#internalPotentialArray = timeStamp(units[4].average, layers, 1)
#eqU = sum(internalPotentialArray[-15:])/len(internalPotentialArray[-15:])

print("surface temperature:", surfaceT[0])
print("Speed of sound:", machSpeedOfSound[1])
print("Speed of sound [m/s]:", calc_real_sos(machSpeedOfSound[1],m,aeps))
print("internal Energy at equillibrium in dimensionless units:", eqU)
print("Average density in same place:",eqn)

# print(shockFront[0])

# linearRegression(shockFront[1],shockFront[0],'ShockFront.png',"[t*]","[x*]")

# todo: animateArray må fikses
#animateArray(units[0], layers, shockFront[0])


#print(plottingMD(files2, names1, tidspunkt))

# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 09:20:16 2021

@author: Tagew
"""
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
import matplotlib.ticker as mtick

#Her lages en klasse for hver av verdiene der den har et navn en snitt matrise og en standardavvik matrise
class Unit:
  def __init__(self, name, average, deviation, standardError, median):
    self.name = name
    self.average = average
    self.deviation = deviation
    self.standardError = standardError
    self.median = median

# FOR AT KODEN SKAL FUNKE MÅ ADRESSENE gis med ''\\'', IKKE ''\''

a = "C:\\Users\\Tagew\\MD1\\Kjøringer\\kjøring1\\2a\\BPROPS.DAT"
b = "C:\\Users\\Tagew\\MD1\\Kjøringer\\kjøring1\\2b\\BPROPS.DAT"
c = "C:\\Users\\Tagew\\MD1\\Kjøringer\\kjøring1\\2c\\BPROPS.DAT"
d = "C:\\Users\\Tagew\\MD1\\Kjøringer\\kjøring1\\2d\\BPROPS.DAT"
e = "C:\\Users\\Tagew\\MD1\\Kjøringer\\kjøring1\\2e\\BPROPS.DAT"
f = "C:\\Users\\Tagew\\MD1\\Kjøringer\\kjøring1\\2f\\BPROPS.DAT"
g = "C:\\Users\\Tagew\\MD1\\Kjøringer\\kjøring1\\2g\\BPROPS.DAT"
h = "C:\\Users\\Tagew\\MD1\\Kjøringer\\kjøring1\\2h\\BPROPS.DAT"

n1 = "C:\\Users\\Tagew\\MD1\\Kjøringer\\LinuxParaleller\\BPROPS1.DAT"
n2 = "C:\\Users\\Tagew\\MD1\\Kjøringer\\LinuxParaleller\\BPROPS2.DAT"
n3 = "C:\\Users\\Tagew\\MD1\\Kjøringer\\LinuxParaleller\\BPROPS3.DAT"
n4 = "C:\\Users\\Tagew\\MD1\\Kjøringer\\LinuxParaleller\\BPROPS4.DAT"
n5 = "C:\\Users\\Tagew\\MD1\\Kjøringer\\LinuxParaleller\\BPROPS5.DAT"
n6 = "C:\\Users\\Tagew\\MD1\\Kjøringer\\LinuxParaleller\\BPROPS6.DAT"
n7 = "C:\\Users\\Tagew\\MD1\\Kjøringer\\LinuxParaleller\\BPROPS7.DAT"
n8 = "C:\\Users\\Tagew\\MD1\\Kjøringer\\LinuxParaleller\\BPROPS8.DAT"
n9 = "C:\\Users\\Tagew\\MD1\\Kjøringer\\LinuxParaleller\\BPROPS9.DAT"
n10 = "C:\\Users\\Tagew\\MD1\\Kjøringer\\LinuxParaleller\\BPROPS10.DAT"

m1 = "C:\\Users\\Tagew\\MD1\\Kjøringer\\LinuxAnemdFiks\\CPROPS1.DAT"
m2 = "C:\\Users\\Tagew\\MD1\\Kjøringer\\LinuxAnemdFiks\\CPROPS2.DAT"
m3 = "C:\\Users\\Tagew\\MD1\\Kjøringer\\LinuxAnemdFiks\\CPROPS3.DAT"
m4 = "C:\\Users\\Tagew\\MD1\\Kjøringer\\LinuxAnemdFiks\\CPROPS4.DAT"
m5 = "C:\\Users\\Tagew\\MD1\\Kjøringer\\LinuxAnemdFiks\\CPROPS5.DAT"
m6 = "C:\\Users\\Tagew\\MD1\\Kjøringer\\LinuxAnemdFiks\\CPROPS6.DAT"
m7 = "C:\\Users\\Tagew\\MD1\\Kjøringer\\LinuxAnemdFiks\\CPROPS7.DAT"
m8 = "C:\\Users\\Tagew\\MD1\\Kjøringer\\LinuxAnemdFiks\\CPROPS8.DAT"
m9 = "C:\\Users\\Tagew\\MD1\\Kjøringer\\LinuxAnemdFiks\\CPROPS9.DAT"
m10 = "C:\\Users\\Tagew\\MD1\\Kjøringer\\LinuxAnemdFiks\\CPROPS10.DAT"


#Filvariablene fylles inn i filer 

filer = [a,b,c,d,e,f,g,h]
files1 = [n1,n2,n3,n4,n5,n6,n7,n8,n9,n10]
files2 = [m1,m2,m3,m4,m5,m6,m7,m8,m9,m10]
#Farge rekkefølgen på plottene

colors = ['b', 'g','r','c' ,'m', 'y','k','b', 'g','r','c' ,'m', 'y','k','b', 'g','r','c' ,'m', 'y','k']

#Navnerekkefølge på dataene du legger inn. (kommer nok til å oppdatere denne delen så man slipper å skrive inn alle variabelene)

names = ['T','n',"p",'x1','up', 'h', 'Jm1','Jm2','ju1','ju2','Tx','Ty','Tz','px','py','pz']
names1 = ['T','n',"p",'x1','up', 'h', 'Jm1','Jm2','ju1','ju2','Tx','Ty','Tz','px','py','pz','vcm', 'jqflow', 'jq' ]



#MINOR FUNCTIONS NEEDED TO RUN THE LARGER FUNCTIONS
def timeStamp(array, cLength, time):
    Y = []
    for i in range(len(cLength)):
                Y.append(array[i][time])
    return Y
def regToList(koeff,inter, xList):
    y = []
    for i in range(len(xList)):
        yi = koeff*xList[i] + inter
        y.append(yi)
    return y
def takeClosest(middleValue,y):
   maxValue = max(y)
   maxPos = y.index(maxValue)
   k = 1000
   for i in range(maxPos,len(y)):
       if abs(y[i] - middleValue) < abs(k - middleValue):
           k = y[i]
           index = i  
   return index
def isNaN(num):
    #sjekker om verdien som leses av er et tall
    return num != num
def columnLength(array):
    #regner ut antall kolonner/lag
    for i in range(len(array)):
        if isNaN(array[i][0]):
            break
    return i
def rowLength(array):
    #lengden på radene, dvs antal tidssteg
    a = len(array[1])-2
    return a
def layerList(array):
    layerList = []
    #henter ut en liste med laglengder
    for i in range(columnLength(array)):
        k = float(array[i][1])
        layerList.append(k)
    return layerList
def noGo(number, variables):
    liste = []
    for i in range(1,variables+1):
        liste.append((i*number)+1)
    return liste
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False
def intervals(times,length):
    liste = [0, length-1]
    for i in range(1,times):
        liste.append(liste[-1]+1)
        liste.append(liste[-1]+length-1)
    return liste
#END OF MINOR FUNCTIONSLIST

#Funskjonen averageData finner antall lage og regner ut snitt og standardavvik for filene
def averageData(filenames):
    #lager tomme 
    values = []
    averages = []
    deviations = []
    errors = []
    medians = []
    # Her hentes ut data fra første filen for å lage en liste over lagene
    data = pd.read_fwf(filenames[0])
    array0 = data.to_numpy()
    layers = layerList(array0)
    for i in filenames:
        #Her leses datanene fra filene av. De blir gjort om til store matriser og legges inn i en liste
        data = pd.read_fwf(i)
        array = data.to_numpy()
        values.append(array)
    for l in range(len(array0)):
            #Her går programmet gjennom fil-matrisene og regner ut snittverdi og standardavvik
            aveRow = []
            sdRow = []
            seRow = []
            medRow = []
            for m in range(2,len(array[l])):
                sample = []
                for k in range(len(filenames)):
                    try:
                        if isNaN(float(values[k][l][m])):
                            #sjekker om det 
                            break
                        else:
                            b = float(values[k][l][m])
                            sample.append(b)
                    except ValueError:
                        print (values[k][l][m])
                        break
                if len(sample) > 0:
                    #regner ut snitt og standard avvik
                    ave = statistics.mean(sample)
                    sd = 3*statistics.stdev(sample)
                    sa = sd/math.sqrt(len(sample))
                    med = statistics.median(sample)
                    aveRow.append(ave)
                    sdRow.append(sd)
                    seRow.append(sa)
                    medRow.append(med)
            if len(sdRow) > 0:
                #etter at snitt og SD er regnet ut for raden legges raden inn i snitt og Sd matrise
                averages.append(aveRow)
                deviations.append(sdRow)
                errors.append(seRow)
                medians.append(medRow)
    #Programmet returnerer en liste med lagavstand, en snittmatrise og en sdmatrise
    del aveRow
    del sdRow
    del seRow
    del medRow
    del array
    del data
    return [layers,averages, deviations, errors, medians]

def processing(filenames, names):
    #Formålet med denne funskjonen er å dele snitt og SD matrisen til en liste klasser
    #der hver klasse er en enhet med sitt respektive navn, snitt og sd matrise
    fileRead = averageData(filenames)
    layerList = fileRead[0]
    averages = fileRead[1]
    deviations = fileRead[2]
    errors = fileRead[3]
    medians = fileRead[4]
    values = []
    cLength = len(layerList)
    interval = intervals(len(names),cLength)
    for i in range(len(names)):
        #Her splittes matrisene opp og settes inn i en klasse
        ave = averages[(interval[2*i]):(interval[(2*i)+1]+1)]
        sd = deviations[(interval[2*i]):(interval[(2*i)+1]+1)]
        se = errors[(interval[2*i]):(interval[(2*i)+1]+1)]
        med = medians[(interval[2*i]):(interval[(2*i)+1]+1)]
        name = names[i]
        p = Unit(name, ave, sd,se,med)
        values.append(p)
    return [layerList,values] #programmet returnerer en liste med lagene og en liste med klassene

def plottingMD(filenames,names, time): 
    #funksjonen tar inn en liste med filer, rekkefølgen på navn og tidspunktet man ønsker å undersøke
    #her hentes dataene ut og plottes
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
                #Her plottes figurene
                plt.errorbar(layers, Y, E, linestyle='None', marker='.',color=colors[i])
                if max(Y) > 0:
                    plt.axis([0,layers[-1]+5,0, (max(Y)+max(Y)*0.1)])
                else:
                    plt.axis([0,layers[-1]+5,0, (max(Y)+0.1)])
                plt.rcParams['axes.spines.right'] = False
                plt.rcParams['axes.spines.top'] = False
                
                
                
                plt.xlabel('x*')
                stringName = 'figures/'+ values[i].name + '10sim.png'
                yName =  values[i].name +'*'
                plt.ylabel(yName)
                plt.savefig(stringName, dpi = 300)
                plt.show()
    return None





def writeFileVersion2(filenames,names):
    process = processing(filenames, names)
    layers = process[0]
    values = process[1]
    #legge til lag og liste i alle 
    for j in range(len(values)):
        averages = values[j].average
        deviations = values[j].deviation
        standardErrors = values[j].standardError
        medians  = values[j].median
        k = 0  
        for i in range(len(values[0].average)):
            if k == 64:
                k = 0
            layer = layers[k]
            N = k + 1
            averages[i].insert(0,layer)
            averages[i].insert(0,N)
            deviations[i].insert(0,layer)
            deviations[i].insert(0,N)
            standardErrors[i].insert(0,layer)
            standardErrors[i].insert(0,N)
            medians[i].insert(0,layer)
            medians[i].insert(0,N)
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
        sdName =  values[m].name + 'standardDeviatons.dat'
        sdNames.append(sdName)
        seName =  values[m].name + 'standardErrors.dat'
        seNames.append(seName)
        medName = values[m].name + 'medians.dat'
        medNames.append(medName)
        head = 'layer \t x* \t' + values[m].name
        np.savetxt(aName,values[m].average , fmt='%.4E',header=head, 
           comments='')
        np.savetxt(sdName,values[m].deviation, fmt='%.4E', header=head, 
           comments='')
        np.savetxt(seName,values[m].standardError, fmt='%.4E', header=head, 
           comments='')
        np.savetxt(medName, values[m].median, fmt='%.4E',header=head, 
           comments='')
    fout=open("Averages.dat","w")
    
    for line in open(aveNames[0]):
        fout.write(line)  
    for num in range(1,len(aveNames)):
        f = open(aveNames[num])
        for line in f:
             fout.write(line)
        f.close()
    fout.close()
    fout=open("StandardDeviations.dat","w")
    for line in open(sdNames[0]):
        fout.write(line)
    for num in range(1,len(aveNames)):
        f = open(sdNames[num])
        for line in f:
             fout.write(line)
        f.close()
    fout.close()
    fout=open("StandardErrors.dat","w")
    for line in open(seNames[0]):
        fout.write(line) 
    for num in range(1,len(aveNames)):
        f = open(seNames[num])
        for line in f:
             fout.write(line)
        f.close()
    fout.close()
    fout=open("Medians.dat","w")
    for line in open(medNames[0]):
        fout.write(line) 
    for num in range(1,len(aveNames)):
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


def findGibbsSurSimp(densities, layers, animation = False, plotting = False):
    #for løkke for relevante tidspunkt
    stop = len(densities[0])-6
    shockPosition = []
    filenames = []
    time = []
    maxYaxis = max(timeStamp(densities,layers, 5))
    for i in range(4,stop):
        time.append(i)
        data = timeStamp(densities,layers,i)
        maxValue = max(data)
        maxPos = data.index(maxValue)
        print(maxPos)
        eqPos = maxPos +20
        yfAve = (sum(data[-10:-1])/len(data[-10:-1]))
        yAveList = []
        for j in range(len(layers)):
            yAveList.append(yfAve)
        cV = (maxValue +data[eqPos])/2
        c = takeClosest(cV,data)
        a = c - 5
        b = c + 5
        aN = 10
        bN = 10
        I1 = simps(data[a-aN:b+bN+1],layers[a-aN:b+bN+1])
        #Making linear regression lines for the front(f) and back(b)
        #kanskje lurere å basere frontlinjen på et snitt av likevektsverdiene
        yb = np.array(data[(a-aN):a+1])
        xb = np.array(layers[(a-aN):a+1]).reshape((-1, 1))
        #yf = np.array(data[b:(b+bN+1)])
        #xf = np.array(layers[b:(b+bN+1)]).reshape((-1, 1))
        
        #Create a model and fit it
        modelBack = LinearRegression()
        modelBack.fit(xb, yb)
        #modelFront = LinearRegression()
        #modelFront.fit(xf,yf)
        
        aBack = modelBack.coef_
        bBack = modelBack.intercept_        
        #aFront = modelFront.coef_
        #bFront = modelFront.intercept_
        yBack = regToList(aBack,bBack,layers)
        #yFront= regToList(aFront,bFront,layers)
       
        def f(l):
            #f = ((l-layers[a-aN])*(yBack[a-aN]+yBack[a])/2+(layers[(b+bN)]-l)*(yfAve)-I1)**2
            f = (0.5*aBack*(l*l-layers[a-aN]**2) + bBack*(l - layers[a-aN]) + yfAve*(layers[b+bN] -l) - I1)
            #B: (((l-layers[a-aN])*(yBack[a-aN]+yBack[a])/2+(layers[(b+bN)]-l)*(yfAve)-I1)**2)
            return f
        
        def m(l):
            #Bjørn:(((l-layers[a-aN])*(yBack[a-aN]+yBack[a])/2+(layers[(b+bN)]-l)*(yfAve)-I1)**2)
            #((0.5*aBack*(l*l-layers[a-aN]**2) + bBack*(l - layers[a-aN]) + yfAve*(layers[b+bN] -l) - I1)**2)
            return (((0.5*aBack*(l*l-layers[a-aN]**2) + bBack*(l - layers[a-aN]) + yfAve*(layers[b+bN] -l) - I1))**2)
        
        guess = layers[maxPos]
        #Root method
        solution = scipy.optimize.root(f,guess)
        #Minimize method
        #solMin = scipy.optimize.minimize_scalar(m, bounds=(0, layers[-1]),method='bounded')
        
        shockPosition.append(solution['x'][0])
        if plotting == True:
                line = solution['x'][0]
                line2 = solution['x'][0]
                plt.rcParams['axes.spines.right'] = False
                plt.rcParams['axes.spines.top'] = False
                plt.xlabel('x*')
                plt.ylabel('n*')
                plt.plot(layers, yBack, 'b')
                plt.plot(layers, yAveList, 'r')
                plt.plot(layers, data, 'g.', label= 'Density')
                plt.axvline(x=line2)
                plt.axvline(x=layers[a-10])
                plt.axvline(x=layers[b+10])
                plt.axis([0,layers[-1],0,(maxYaxis+maxYaxis*0.25)])
                stringName = 'figures/gibbsEqSur/' + str(i) + 'timeStep.png'
                filenames.append(stringName)
                plt.savefig(stringName, dpi = 300)
                plt.show()
                l = solution['x'][0]
                
    if animation == True and plotting == True:
        with imageio.get_writer('figures/gibbsEqSur/movie.gif', mode='I', duration=1) as writer:
            for filename in filenames:
                image = imageio.imread(filename)
                writer.append_data(image)
        
    return [shockPosition, time]


def linearRegression(xList,yList, filename):
    yb = np.array(yList)
    xb = np.array(xList).reshape((-1, 1))
    modelBack = LinearRegression()
    modelBack.fit(xb, yb)
    aBack = modelBack.coef_
    bBack = modelBack.intercept_
    yBack = regToList(aBack,bBack,xList)
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False
    plt.axis([0,xList[-1]+ xList[-1]*0.1,0,(max(yList)+max(yList)*0.25)])
    plt.plot(xList,yList , 'b.')
    plt.plot(xList, yBack)
    plt.savefig(filename, dpi = 300)
    plt.show()
    return None

def animateArray(classType, layers, shockPosition):
    array = classType.average
    se = classType.standardError
    maxYaxis = max(timeStamp(array,layers, 4))
    
    filenames = []
    for i in range(len(array[0])):
        yList = timeStamp(array, layers, i)
        eList = timeStamp(se, layers,i )
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
                        
                        
        plt.axis([0,layers[-1]+0.1*layers[-1],0,(maxYaxis+maxYaxis*0.25)])
        stringName = 'figures/animate/' + classType.name + str(i) + 'timeStep.png'
        filenames.append(stringName)
        plt.savefig(stringName, dpi = 300)
        plt.show()
    gifName = 'figures/animate/'+ classType.name +'movie.gif'
    with imageio.get_writer(gifName, mode='I', duration=1) as writer:
            for filename in filenames:
                image = imageio.imread(filename)
                writer.append_data(image)
                
    for m in filenames:
        os.remove(m)
    return None
        

#SETT HVILKET TIDSPUNKT DU VIL KJØRE            
tidspunkt = 10           
#Skriver snitt, SD, SE og Middelverdi til DAT.filer

#writeFileVersion1(filer, names)
#writeFileVersion2(filer,names)


#ha et sett med klasser
#writeFileVersion2(files2,names1)
process = processing(files2, names1)
layers = process[0]
units = process[1]
density = units[1].average


shockFront = findGibbsSurSimp(density,layers)
print(shockFront[0])

#linearRegression(shockFront[1],shockFront[0],'ShockFront.png')

#animateArray må fikses
animateArray(units[16], layers, shockFront[0])


#print(plottingMD(files2, names1, tidspunkt))


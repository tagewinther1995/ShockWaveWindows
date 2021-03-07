# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 14:01:16 2021

@author: Tagew
"""


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

def writeFileVersion1(filenames,names):
    b = averageData(filenames)
    layers = b[0]
    averages = b[1]
    deviations = b[2]
    standardErrors = b[3]
    medians = b[4]
    k = 0
    stringList = ["layer", "x*",names[0]]
    for i in range(len(averages)):
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
    np.savetxt('averages.dat',averages , fmt='%.4E')
    np.savetxt('standardDeviatons.dat',deviations , fmt='%.4E')
    np.savetxt('standardErrors.dat',standardErrors , fmt='%.4E')
    np.savetxt('medians.dat', medians, fmt='%.4E')
    return None
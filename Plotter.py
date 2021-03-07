# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 09:26:28 2021

@author: Tagew
"""

import xlrd 
import numpy as np
import matplotlib.pyplot as plt


def layer_list(nLayers): #Creates a list of all the layers in the Md Simulation
    loc = (r"C:\Users\Tagew\MD1\ExcelData\temps.xlsx")
    wb = xlrd.open_workbook(loc)
    sheet = wb.sheet_by_index(0)
    layer_list = []
    for i in range(nLayers):
        layer_list.append(sheet.cell_value(i+1, 0))
    return layer_list


def layer_list_temp(nLayers, time): #Creates a list of the values for the system at a given time step Simulation
    loc = (r"C:\Users\Tagew\MD1\ExcelData\temps.xlsx")
    wb = xlrd.open_workbook(loc)
    sheet = wb.sheet_by_index(0)
    layer_list = []
    for i in range(nLayers):
        layer_list.append(float(sheet.cell_value(i+1, time)))
    return layer_list


def layer_list_dens(nLayers, time): #Creates a list of the values for the system at a given time step Simulation
    loc = (r"C:\Users\Tagew\MD1\ExcelData\dens.xlsx")
    wb = xlrd.open_workbook(loc)
    sheet = wb.sheet_by_index(0)
    layer_list = []
    for i in range(nLayers):
        layer_list.append(float(sheet.cell_value(i+1, time)))
    return layer_list

def layer_list_mflux(nLayers, time): #Creates a list of the values for the system at a given time step Simulation
    loc = (r"C:\Users\Tagew\MD1\ExcelData\mflux.xlsx")
    wb = xlrd.open_workbook(loc)
    sheet = wb.sheet_by_index(0)
    layer_list = []
    for i in range(nLayers):
        layer_list.append(float(sheet.cell_value(i+1, time)))
    return layer_list

def layer_list_press(nLayers, time): #Creates a list of the values for the system at a given time step Simulation
    loc = (r"C:\Users\Tagew\MD1\ExcelData\press.xlsx")
    wb = xlrd.open_workbook(loc)
    sheet = wb.sheet_by_index(0)
    layer_list = []
    for i in range(nLayers):
        layer_list.append(float(sheet.cell_value(i+1, time)))
    return layer_list

time = 11
N = 64

def plot_list(N, time): #Establishes a plotting list 
    layers = layer_list(N)
    temps = layer_list_temp(N, time)
    dens = layer_list_dens(N, time)
    mflux = layer_list_mflux(N, time)
    press = layer_list_press(N, time)
    return [layers, temps, dens, mflux, press]



tot_list = plot_list(N,time)

plt.plot(tot_list[0], tot_list[1], '.', label= 'Temperature')
plt.xlabel('Layer Number')
plt.ylabel('Temperature[T*]')
plt.savefig('Temperature.png', dpi = 300)
plt.show()

plt.plot(tot_list[0], tot_list[2], 'r--', label= 'Density')
plt.xlabel('Layer Number')
plt.ylabel('Density[n*]')
plt.savefig('Density.png', dpi = 300)
plt.show()

plt.plot(tot_list[0], tot_list[3], 'g--', label= 'Mass flux')
plt.xlabel('Layer Number')
plt.ylabel('Massflux[Jm*]')
plt.savefig('mflux.png', dpi = 300)
plt.show()

plt.plot(tot_list[0], tot_list[4], 'c--', label= 'pressure')
plt.xlabel('Layer Number')
plt.ylabel('Pressure[p*]')
plt.savefig('press.png', dpi = 300)
plt.show()

plt.plot(tot_list[0], tot_list[1], 'b--', label= 'Temperature')
plt.plot(tot_list[0], tot_list[2], 'r--', label= 'Density')
plt.plot(tot_list[0], tot_list[3], 'g--', label= 'Mass flux')
plt.plot(tot_list[0], tot_list[4], 'c--', label= 'pressure')
plt.xlabel('Layer Number')
plt.ylabel('Dimless unit')
plt.legend()
plt.savefig('Total.png', dpi = 300)
plt.show()





# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 12:40:34 2021

@author: Tagew
"""

import scipy.optimize

def coffein_extraction (unknowns):
    """
    Input : List of unknown variables [n2 , n4 , n5]
    """
    # Unpacking unknowns :
    m3 , m4, x3V = unknowns
    # Defining known variables
    # Stream 1:
    x1O = 0 # [-]
    x1V = 100/100.07 # [-]
    x1K = 0.07/100.07 # [-]
    m1 = 100.07 # [g]
    # Stream 2:
    x2O = 1.0 #[-]
    x2V = 0 #[-]
    x2K = 0 #[-]
    m2 = 50 #[g]
    # Stream 3:
    x3O = 0.00 # [-]
    #x3V = ? # [-]
    x3K = 1-x3V # [-]
    #m3 = ? # [g]
    4
    # Stream 4:
    x4O = 0.999 # [-]
    x4V = 0.0 # [-]
    x4K = 0.001 # [-]
    #m4 ? # [g]
    # Defining mass balance equations as f1 , f2 , f3
    f1 = m1 * x1O + m2 * x2O - m3 * x3O - m4 * x4O #Mass balance CH2Cl2
    f2 = m1 * x1V + m2 * x2V - m3 * x3V - m4 * x4V #Mass balance water
    f3 = m1 + m2 - m3 - m4 #Mass balance total
    return [f1 , f2 , f3]
# Set initial guesses for unknown values :
m3_guess = 100
m4_guess = 100
x3V_guess = 0.9
m_guess = [m3_guess , m4_guess , x3V_guess]
solution = scipy.optimize.root (coffein_extraction, m_guess)
print('m3 = ' , '%.2f' %solution['x'][0])
print('m4 = ' , '%.2f' %solution['x'][1])
print('x3V = ' ,'%.4f' %solution['x'][2])
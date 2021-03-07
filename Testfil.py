# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 13:26:15 2021

@author: Tagew
"""
import numpy as np
class unit:
  def __init__(self, name, matrix, deviations):
    self.name = name
    self.matrix = matrix
    self.deviations = deviations

a = np.array([1,2,3])
b = np.array([1,2,3])
c = [a,b]
q = [1,2,3,4,5]
a = np.array([[1, 2, 3], ['s', '5', ' '], [7, 8, 9]])
mat = np.array([a])
with open('outfile.dat','wb') as f:
    for line in mat:
        try:
            np.savetxt(f, line, fmt='%.2f')
        except TypeError:
            np.savetxt(f, line, fmt='%s')
print(q[0])
print(q[1])
q.insert(0,88)
print(q[0])
print(q[1])

streng = 'haha'

p = unit('roar',[2131],[2313])
p.name = "geir"
print(p.name)
averages = np.array([[1,2,3,4],[1,2,3,4]])

np.savetxt('test.dat',averages , fmt='%.4E',header="ID \t AMOUNT", 
           comments='')

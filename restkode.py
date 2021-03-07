# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 14:31:33 2021

@author: Tagew
"""

space = ' '
    b = 0
    for i in range(len(names)):
        head = ['layer','x*',names[i]]
        for j in range(len(averages[0])-3):
            head.append(space)
        averages.insert((i*len(layers))+b,head)
        b = b+1
             
    mat = np.array([averages])
    with open('outfile.dat','wb') as f:
        for line in mat:
            try:
                np.savetxt(f, line, fmt='%.4E')
            except TypeError:
                np.savetxt(f, line, fmt='%s')
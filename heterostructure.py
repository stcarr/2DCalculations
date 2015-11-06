# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 18:57:44 2015

@author: Stephen
"""

class Heterostructure:
    
    def __init__(self,sheets,angles,heights):
        self.sheets = sheets
        self.angles = angles
        self.heights = heights
        self.setIndex()
        
    def setIndex(self):
        
        k_total = 0
        s_total = 0
        
        for s in self.sheets:
            k_total += s.max_index
            s_total += 1
        
        self.max_index = k_total
        self.max_sheets = s_total
        
    def posAtomIndex(self, k):
        if k >= self.max_index:
            raise RuntimeError('Hamiltonian Index ' + str(k) + ' is out of bounds in Heterostructure!')
        
        findSheet = True
        current_sheet = 0
        current_index = 0
        
        while findSheet:
            if k < (current_index + self.sheets[current_sheet].max_index):
                sheet_index = current_sheet
                findSheet = False
            else:
                current_index += self.sheets[current_sheet].max_index
                current_sheet += 1
                
        local_pos = self.sheets[sheet_index].posAtomIndex(k - current_index)
        local_x = local_pos[0]
        local_y = local_pos[1]
        local_z = local_pos[2]
        theta = self.angles[sheet_index]
        x = local_x*math.cos(theta) - local_y*math.sin(theta)
        y = local_y*math.cos(theta) + local_x*math.sin(theta)
        z = local_z + self.heights[sheet_index]
        
        return (x,y,z)
        
    def posAtomGrid(self,grid_index,s):
        
        local_pos = self.sheets[s].posAtomGrid(grid_index)
        local_x = local_pos[0]
        local_y = local_pos[1]
        local_z = local_pos[2]
        theta = self.angles[s]
        x = local_x*math.cos(theta) - local_y*math.sin(theta)
        y = local_y*math.cos(theta) + local_x*math.sin(theta)
        z = local_z + self.heights[s]
        
        return (x,y,z)
        
    def find_nearest(self,pos,s):
        theta = self.angles[s]
        
        x = pos[0]
        y = pos[1]
        z = pos[2]
        
        x_new = x*math.cos(theta) + y*math.sin(theta)
        y_new = y*math.cos(theta) - x*math.sin(theta)
    
        i = min(max(int(math.floor(x_new)),self.sheets[s].min_shape[0]),self.sheets[s].max_shape[0]) - self.sheets[s].min_shape[0] 
        j = min(max(int(math.floor(y_new)),self.sheets[s].min_shape[1]),self.sheets[s].max_shape[1]) - self.sheets[s].min_shape[1]
        
        return (i,j,0)
        
        

from sheet import Sheet

import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


if __name__ == '__main__':
    N = 5
    s1 = Sheet([[1,0,0],[0,1,0],[0,0,1]],['C','C'],[[0,0,0],[0.5,0.5,0]],[0],[-N,-N],[N,N],'graphene')
    s2 = Sheet([[1,0,0],[0,1,0],[0,0,1]],['C','C'],[[0,0,0],[0.5,0.5,0]],[0],[-N,-N],[N,N],'graphene')
    h = Heterostructure((s1,s2),[0,math.pi/60],[0,6])

    k_max = h.max_index 
    pos_array_x = []
    pos_array_y = []
    pos_array_z = []
    for k in xrange(k_max):
        pos_array_x.append(h.posAtomIndex(k)[0])
        pos_array_y.append(h.posAtomIndex(k)[1])
        pos_array_z.append(h.posAtomIndex(k)[2])
        
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(pos_array_x,pos_array_y,pos_array_z)
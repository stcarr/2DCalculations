# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 17:10:21 2015

@author: Stephen
"""

import time
import numpy as np
import matplotlib.pyplot as plt

from scipy.sparse import lil_matrix
from intralayers import Intralayer

class Sheet:
    

# Defines one sheet's geometery for a 2d heterostructure simulation
# Takes in variables:
#
#   unitcell = [a1,a2,a3], a vector of vectors which defines the unit cell in real space
#   atom_types = ['C','N',...] vector of chacters of atomic symbol
#   atom_pos = [pos1,pos2,...] vector of atomic positions in unit cell
#   shape = input function for the shape of the sheet (not yet implemented)
#   min_shape = [min_i,min_j] sets up lower cutoff for our grid search over shape
#   max_shape = [max_i,max_j] sets up upper cutoff for our grid search over shape
#
#  ///
    
    def __init__(self, unitcell, atom_types, atom_pos, shape, min_shape, max_shape, mat):
        self.a1 = unitcell[0]
        self.a2 = unitcell[1]
        self.a3 = unitcell[2]
        self.atom_types = atom_types
        self.atom_pos = atom_pos        
        self.shape = shape
        self.mat = mat
        
        # max shape should be a tuple (i,j) = (max rows, max cols)
        self.max_shape = max_shape
        self.min_shape = min_shape
        self.setIndex()
    
    # setup for the sheet's geometery and indexing
    def setIndex(self):
        
        self.grid_array = [[[-1 for s in xrange(len(self.atom_types))] for j in xrange(self.min_shape[1],self.max_shape[1])] for i in xrange(self.min_shape[0],self.max_shape[0])]       
        k = 0
        
        for i in xrange(self.max_shape[0] - self.min_shape[0]): 
            for j in xrange(self.max_shape[1] - self.min_shape[1]):
                for l in xrange(len(self.atom_types)):
                    if self.checkShape(self.posAtomGrid([i,j,l])) == True:
                        self.grid_array[i][j][l] = k
                        k += 1
        
        self.max_index = k
        
        self.index_array = [[0,0,0] for m in xrange(self.max_index)]
        
        for i in xrange(self.max_shape[0] - self.min_shape[0]): 
            for j in xrange(self.max_shape[1] - self.min_shape[1]):
                for l in xrange(len(self.atom_types)):
                    k_here = self.grid_array[i][j][l]
                    if k_here != -1:
                        self.index_array[k_here] = [i,j,l]
                                
    # check to see if a position is in the sheet as defined by shape                  
    def checkShape(self, pos):
        #currently just makes the sheet a circle of radius = max_shape[0]        
        
        x = pos[0]
        y = pos[1]
        z = pos[2]
        
        if x**2 + y**2 < self.max_shape[0]**2:
            return True
        else:
            return False
   
     
    # position of atom at index k                
    def posAtomIndex(self, k):
        return self.posAtomGrid(self.indexToGrid(k))
    
    # position of atom s at unit cell (i,j)
    def posAtomGrid(self, grid_index):
        i = grid_index[0] + self.min_shape[0]
        j = grid_index[1] + self.min_shape[1]
        l = grid_index[2]        
        
        x = float(i)*self.a1[0]+float(j)*self.a2[0] + self.atom_pos[l][0]
        y = float(i)*self.a1[1]+float(j)*self.a2[1] + self.atom_pos[l][1]
        z = self.atom_pos[l][2]
        
        return (x,y,z)
    
    # convert index k to grid position (i,j,s)
    def indexToGrid(self, k):
    
        if k >= self.max_index:
            raise RuntimeError('Hamiltonian Index ' + str(k) + ' is out of bounds in Sheet!')
                
        
        i = self.index_array[k][0]
        j = self.index_array[k][1]
        s = self.index_array[k][2]
        
        return [i,j,s]
      
    # convert grid position (i,j,s) to index k      
    def gridToIndex(self, grid_index):
        i = grid_index[0]
        j = grid_index[1]
        l = grid_index[2]  
        
        k = self.grid_array[i][j][l]
        return k
        
    def intraHamiltonian(self,type):
        print 'Calculating intralyer H...'
        H = lil_matrix((self.max_index,self.max_index))
        
        for k in xrange(self.max_index):
            
            grid_here = self.indexToGrid(k)
            ilay = Intralayer()
            energies = ilay.interaction(self.mat,grid_here)
            
            for l in xrange(0,len(energies)):
                i = energies[l][0]
                j = energies[l][1]
                s = energies[l][2]
                t = energies[l][3]
                
                if i > 0 and i < self.max_shape[0] - self.min_shape[0]:
                    if j > 0 and j < self.max_shape[1] - self.min_shape[1]:
                        if s > 0 and s < self.atom_types:
                            new_index = self.gridToIndex([i,j,s])
                            
                            if new_index != -1:
                                H[k,new_index] = t
            
        return H.tocsr()

import matplotlib.pyplot as plt

if __name__ == '__main__':
    N_array  = [5,10,20,30,40,50]
    t_array  = [0, 0, 0, 0, 0, 0]
    i = 0
    for N in N_array:
        tic = time.time()
        s = Sheet([[1,0,0],[0,1,0],[0,0,1]],['C','C'],[[0,0,0],[0.5,0.5,0]],[0],[-N,-N],[N,N],'graphene')
        
        H = s.intraHamiltonian(0)
        print H.getnnz
        
        t_array[i] = time.time() - tic
        i += 1
        
    p = np.poly1d(np.polyfit(N_array,t_array,2))
    plt.plot(N_array,t_array,'b')
    plt.plot(xrange(0,1000),p(xrange(0,1000)),'r')
    
    '''
    k_max = s.max_index
    pos_array_x = []
    pos_array_y = []
    for k in xrange(k_max):
        pos_array_x.append(s.posAtomIndex(k)[0])
        pos_array_y.append(s.posAtomIndex(k)[1])
        
    plt.scatter(pos_array_x,pos_array_y)
    '''
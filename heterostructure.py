# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 18:57:44 2015

@author: Stephen
"""

from interlayers import Interlayer
from locality import Locality
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg
import numpy as np
import matplotlib.pyplot as plt
import time

class Heterostructure:
    
    #heterostructure is defined as a collection of sheets, angles, and heights
    def __init__(self,sheets,angles,heights):
        self.sheets = sheets
        self.angles = angles
        self.heights = heights
        self.setIndex()
        shift = []        
        for i in sheets:
            shift.append((0,0,0))
            
        self.shift = shift
    
    #counts the global index (for making a complete H matrix)
    def setIndex(self):
        
        k_total = 0
        s_total = 0
        
        for s in self.sheets:
            k_total += s.max_index
            s_total += 1
        
        self.max_index = k_total
        self.max_sheets = s_total
        
    def setShift(self,sheet,b):
        self.shift[sheet] = b
    
    #returns an atom's poisition given its (global) index k
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
        
        s = sheet_index
        local_pos = self.sheets[s].posAtomIndex(k - current_index)
        local_x = local_pos[0]
        local_y = local_pos[1]
        local_z = local_pos[2]
        theta = self.angles[s]
        
        x = local_x*math.cos(theta) - local_y*math.sin(theta) + self.shift[s][0]*self.sheets[s].a1[0] + self.shift[s][1]*self.sheets[s].a2[0] + self.shift[s][2]*self.sheets[s].a3[0]
        y = local_y*math.cos(theta) + local_x*math.sin(theta) + self.shift[s][0]*self.sheets[s].a1[1] + self.shift[s][1]*self.sheets[s].a2[1] + self.shift[s][2]*self.sheets[s].a3[1]
        z = local_z + self.heights[s] + self.shift[s][0]*self.sheets[s].a1[2] + self.shift[s][1]*self.sheets[s].a2[2] + self.shift[s][2]*self.sheets[s].a3[2]
                
        return (x,y,z)
    
    #gives an atom's position given its grid postion on a sheet
    def posAtomGrid(self,grid_index,s):
        
        local_pos = self.sheets[s].posAtomGrid(grid_index)
        local_x = local_pos[0]
        local_y = local_pos[1]
        local_z = local_pos[2]
        theta = self.angles[s]
        
        x = local_x*math.cos(theta) - local_y*math.sin(theta) + self.shift[s][0]*self.sheets[s].a1[0] + self.shift[s][1]*self.sheets[s].a2[0] + self.shift[s][2]*self.sheets[s].a3[0]
        y = local_y*math.cos(theta) + local_x*math.sin(theta) + self.shift[s][0]*self.sheets[s].a1[1] + self.shift[s][1]*self.sheets[s].a2[1] + self.shift[s][2]*self.sheets[s].a3[1]
        z = local_z + self.heights[s] + self.shift[s][0]*self.sheets[s].a1[2] + self.shift[s][1]*self.sheets[s].a2[2] + self.shift[s][2]*self.sheets[s].a3[2]
        
        return (x,y,z)
    
    #finds the grid pos of the closest atom on sheet index s (ignores z positions)
    def findNearest(self,pos,s):
        theta = self.angles[s]
        
        x = pos[0]
        y = pos[1]
        z = pos[2]
        
        x_new = x*math.cos(theta) + y*math.sin(theta) - self.shift[s][0]*self.sheets[s].a1[0] + self.shift[s][1]*self.sheets[s].a2[0] + self.shift[s][2]*self.sheets[s].a3[0]
        y_new = y*math.cos(theta) - x*math.sin(theta) - self.shift[s][0]*self.sheets[s].a1[1] + self.shift[s][1]*self.sheets[s].a2[1] + self.shift[s][2]*self.sheets[s].a3[1]
    
        i = min(max(int(math.floor(x_new)),self.sheets[s].min_shape[0]),self.sheets[s].max_shape[0]) - self.sheets[s].min_shape[0] 
        j = min(max(int(math.floor(y_new)),self.sheets[s].min_shape[1]),self.sheets[s].max_shape[1]) - self.sheets[s].min_shape[1]
        
        return (i,j,0)
      
    #returns the interlayer coupling hamiltonian (i.e. H_ij for i != j)
    def interHamiltonian(self,s1,s2):
        sheet1 = self.sheets[s1]
        sheet2 = self.sheets[s2]
        
        searchsize = 4
        
        print 'Calculating interlayer H...'
        H = sparse.lil_matrix((sheet1.max_index,sheet2.max_index))
        
        for k in xrange(sheet1.max_index):
            
            grid_here = sheet1.indexToGrid(k)
            ih = grid_here[0]
            jh = grid_here[1]
            sh = grid_here[2]            
            
            pos_here = self.posAtomGrid(grid_here,s1)
            
            
            # particular for honeycomb lattice
            # finds neighbors to determine lattice angles in graphene bilayer
            if sh == 0:
                neighbors = [[ih,jh,1],[ih,jh-1,1],[ih-1,jh,1]]
            else:
                neighbors = [[ih,jh,0],[ih+1,jh,0],[ih,jh+1,0]]
                
            pos_bond = (0,0,0)
            for n in neighbors:
                if sheet1.gridToIndex(n) != -1:
                     pos_bond = self.posAtomGrid(n,s1)
                     
            a1 = math.sqrt((pos_here[0] - pos_bond[0])**2 + (pos_here[1] - pos_bond[1])**2)

            
            new_grid = self.findNearest(pos_here,s2) 
            
            i0 = new_grid[0]
            j0 = new_grid[1]
            s0 = new_grid[2]
            
            ilay = Interlayer()
            
            for i in xrange(max(0,i0 - searchsize), min(sheet2.max_shape[0] - sheet2.min_shape[0],i0 + searchsize)):
                for j in xrange(max(0,j0 - searchsize), min(sheet2.max_shape[1] - sheet2.min_shape[1],j0 + searchsize)):
                    for s in xrange(len(sheet2.atom_types)):
                        
                        new_grid = [i,j,s]
                        new_index = sheet2.gridToIndex(new_grid)
                        
                        if new_index != -1:
                            new_pos = self.posAtomGrid(new_grid,s2)
                            
                            #still particular to honeycomb lattice
                            if s == 0:
                                neighbors = [[i,j,1],[i,j-1,1],[i-1,j,1]]
                            else:
                                neighbors = [[i,j,0],[i+1,j,0],[i,j+1,0]]
                                
                            new_pos_bond = (0,0,0)
                            for n in neighbors:
                                if sheet2.gridToIndex(n) != -1:
                                     new_pos_bond = self.posAtomGrid(n,s2)
                            
                            a2 = math.sqrt((new_pos[0] - new_pos_bond[0])**2 + (new_pos[1] - new_pos_bond[1])**2)    
                                 
                            r = math.sqrt((new_pos[0] - pos_here[0])**2 + (new_pos[1] - pos_here[1])**2)
                            r12 = math.sqrt((pos_bond[0] - new_pos[0])**2 + (pos_bond[1] - new_pos[1])**2)
                            r21 = math.sqrt((new_pos_bond[0] - pos_here[0])**2 + (new_pos_bond[1] - pos_here[1])**2)
                            
                            energy = ilay.interaction(sheet1.mat,sheet2.mat,grid_here,new_grid,r,a1,a2,r12,r21)
                            H[k,new_index] = energy
            
        return H.tocsr()
        
    def totalHamiltonian(self):
        # creates an array of sparse matrices, "rows", which can be used to 
        # generate a total hamiltonian by using the scipy sparse.bmat method.        
        
        intraH = []
        interH = []
        rows = []
        smax = len(self.sheets)
        for s in xrange(smax):
            intraH.append(self.sheets[s].intraHamiltonian(0))
            if s != smax - 1:
                interH.append(self.interHamiltonian(s,s+1))
            
            if s == 0:
                row = []
                row.append(intraH[s])
                row.append(interH[s])
                for i in xrange(s+2, smax):
                    row.append(None)
            elif s == smax - 1:
                row = []
                for i in xrange(0, s - 1):
                    row.append(None)
                row.append(sparse.csr_matrix.transpose(interH[s-1]))
                row.append(intraH[s])
            else:
                row = []
                for i in xrange(0,s-1):
                    row.append(None)
                row.append(sparse.csr_matrix.transpose(interH[s-1]))
                row.append(intraH[s])
                row.append(interH[s])
                for j in xrange(s+2,smax):
                    row.append(None)
            
            rows.append(row)
        return rows
        
    def localitySweep(self,samples):
        # implements locality over all sheets in a sample x sample grid of
        # each sheet's unit cell
        
        L = Locality()
        rows = self.totalHamiltonian()
        maxsheet = len(self.sheets) - 1
        
        self.recursiveSweep(samples,maxsheet,rows,L)
        
        val_arr = L.getVals()
        vals = []
        for i in xrange(len(val_arr)):
            for j in xrange(len(val_arr[i])):
                vals.append(val_arr[i][j])
                
        plt.hist(vals,bins=100)
            
        
            
    def recursiveSweep(self,samples,currentsheet,rows,L):
        # recursive segment of localitySweep
        # goes over all possible shifts across all sheets (except for the bottom sheet)
        # and updates the hamiltoninian in an efficient manner before solving
        # the eigenvalue problem. All results stored in the L locality object
        
        if currentsheet != 0:
            for i in np.linspace(0,1,samples,endpoint=False):
                for j in np.linspace(0,1,samples,endpoint=False):
                    self.setShift(currentsheet,[i,j,0])
                    rows_new = rows
                    if currentsheet != len(self.sheets) - 1:
                        rows_new = self.hUpdate(rows,currentsheet)
                    self.recursiveSweep(samples,currentsheet - 1,rows_new,L)
        else:
            
            rows_fin = self.hUpdate(rows,0)
            H = sparse.bmat(rows_fin).tocsr() 
            vals, vecs = linalg.eigs(H,sigma=2, k=H.shape[0] - 2,which='LM',tol=0.001)
            print "Eigen problem solved"
            L.addVals(vals)
            L.addVecs(vecs)
            
    def hUpdate(self,rows,s):
        # update the block [s, s+1]
        #       (H is symmetric so only need to do one calculation for each pair as [k,k`] = [k',k])
        # go over all non-zero matrix elements in these blocks and update them using interlayers.py
        
        inter = self.interHamiltonian(s,s+1)
        rows_new = rows
        rows_new[s][s+1] = inter
        rows_new[s+1][s] = sparse.csr_matrix.transpose(inter)       
        return rows_new
        
        

from sheet import Sheet

import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


if __name__ == '__main__':
    N = 15

    a = 2.46
    a_1 = [a,0,0]
    a_2 = [a/2.0,a*math.sqrt(3)/2,0]
    sites = [[0.0,0.0,0.0],[a/2.0,a/(2.0*math.sqrt(3)),0.0]]
    
    
    s1 = Sheet([a_1,a_2,[0,0,1]],['C','C'],sites,[0],[-N,-N],[N,N],'graphene')
    s2 = Sheet([a_1,a_2,[0,0,1]],['C','C'],sites,[0],[-N,-N],[N,N],'graphene')
    s3 = Sheet([a_1,a_2,[0,0,1]],['C','C'],sites,[0],[-N,-N],[N,N],'graphene')
    h = Heterostructure((s1,s2,s3),[0,math.pi/90,math.pi/45],[0,6,12])
    H = h.localitySweep(2)
    '''
    val_arr = set()
    vec_arr = set()
    for i in np.linspace(0,1,4,endpoint=False):
        for j in np.linspace(0,1,4,endpoint=False):
            h.setShift([[0,0,0],[i,j,0]])
            H12 = h.interHamiltonian(0,1)
            H21 = h.interHamiltonian(1,0)
            H1 = h.sheets[0].intraHamiltonian(0)
            H2 = h.sheets[1].intraHamiltonian(0)
            H = sparse.bmat([[H1,H12],[H21,H2]]).tocsr()

    
            start = time.time()
            vals, vecs = linalg.eigs(H,sigma=2, k=H.shape[0] - 2,which='LM',tol=0.001)
            print time.time() - start
            #val_arr.add(vals)
            #vec_arr.add(vecs)
    '''
    #vals, vecs = linalg.eigs(H,sigma=2, k=H.shape[0] - 2,which='LM',tol=0.001)
    #plt.hist(vals,bins=25)
    #plt.draw()

    #print vals
    '''
    #PLOTS THE ATOMS IN 3D SPACE
    
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
    '''
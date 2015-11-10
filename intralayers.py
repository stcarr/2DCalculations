# -*- coding: utf-8 -*-
"""
Created on Fri Nov 06 15:16:24 2015

@author: Stephen
"""

class Intralayer:
        
        def __init__(self):
            null = 0
            
        def interaction(self,mat,grid):
            if mat == 'graphene':
                return self.int_graphene(grid)
            
        def int_graphene(self,grid):
            
            #!! need to look up which t's correspond to which pair geometery!!
            
            t1 = -2.892
            t2 =  0.243
            t3 = -0.266
            t4 =  0.024
            t5 =  0.052
            t6 = -0.021
            t7 = -0.015
            t8 = -0.021
            
            out_array = [(0,0,0,0) for s in xrange(0,21)]
            k = 0            
            
            for n1 in self.n_graphene(grid):
                
                i1 = n1[0]
                j1 = n1[1]
                s1 = n1[2]
                
                out_array[k] = [i1,j1,s1,t1]
                k += 1
                
                for n2 in self.n_graphene(n1):
                    
                    if n2 != grid:
                
                        i2 = n2[0]
                        j2 = n2[1]
                        s2 = n2[2]
                        
                        out_array[k] = [i2,j2,s2,t2]
                        k += 1
                        
                        for n3 in self.n_graphene(n2):
                            
                            if n3 != n1:
                                                
                                i3 = n3[0]
                                j3 = n3[1]
                                s3 = n3[2]
                                
                                out_array[k] = [i3,j3,s3,t3]
                                k += 1
            return out_array     
                
        def n_graphene(self,grid):
            i = grid[0]
            j = grid[1]
            s = grid[2]
            
            if s == 0:
                return [[i,j,1],[i,j-1,1],[i+1,j,1]]
            if s == 1:
                return [[i,j,0],[i-1,j,0],[i,j+1,0]]
                
if __name__ == '__main__':
    ilay = Intralayer()
    energies = ilay.interaction('graphene',[15,15,0])
    print len(energies)
    print energies[17][0]
    print energies[17][1]
    print energies[17][2]
    print energies[17][3]
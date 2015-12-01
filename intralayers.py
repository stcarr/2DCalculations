# -*- coding: utf-8 -*-
"""
Created on Fri Nov 06 15:16:24 2015

@author: Stephen
"""

class Intralayer:
        
        # Initialising Intralayer does not change any of its properties
        # this is used as a class which returns intralayer hamiltonian energies
        # and so it is used as a sort of library look-up by sheet.py
        def __init__(self):
            null = 0
            
        def interaction(self,mat,grid):
            if mat == 'graphene':
                return self.int_graphene(grid)
            
        def int_graphene(self,grid):
            
            #!! need to look up which t's correspond to which pair geometery!!
                        
            i = grid[0]
            j = grid[1]
            s = grid[2]
            
            t1 = -2.892
            t2 =  0.243
            t3 = -0.266
            t4 =  0.024
            t5 =  0.052
            t6 = -0.021
            t7 = -0.015
            t8 = -0.021
            
            outarray = [(0,0,0,0) for y in xrange(0,39)]
                
            
            if s == 0:
                #t1 hopping
                outarray[0]  = (i  ,j  ,1,t1)
                outarray[1]  = (i-1,j  ,1,t1)
                outarray[2]  = (i  ,j-1,1,t1)
        
                #t2 hopping
                outarray[3]  = (i  ,j+1,0,t2)
                outarray[4]  = (i+1,j  ,0,t2)
                outarray[5]  = (i+1,j-1,0,t2)
                outarray[6]  = (i  ,j-1,0,t2)
                outarray[7]  = (i-1,j  ,0,t2)
                outarray[8]  = (i-1,j+1,0,t2)
                
                #t3 hopping
                outarray[9]  = (i-1,j-1,1,t3)
                outarray[10] = (i+1,j-1,1,t3)
                outarray[11] = (i-1,j+1,1,t3)
                
                #t4 hopping
                outarray[12] = (i  ,j+1,1,t4)
                outarray[13] = (i+1,j  ,1,t4)
                outarray[14] = (i+1,j-2,1,t4)
                outarray[15] = (i  ,j-2,1,t4)
                outarray[16] = (i-2,j  ,1,t4)
                outarray[17] = (i-2,j+1,1,t4)
                
                #t5 hopping
                outarray[18] = (i-1,j+2,0,t5)
                outarray[19] = (i+1,j+1,0,t5)
                outarray[20] = (i+2,j-1,0,t5)
                outarray[21] = (i+1,j-2,0,t5)
                outarray[22] = (i-1,j-1,0,t5)
                outarray[23] = (i-2,j+1,0,t5)
                
                #t6 hopping
                outarray[24] = (i  ,j+2,0,t6)
                outarray[25] = (i+2,j  ,0,t6)
                outarray[26] = (i+2,j-2,0,t6)
                outarray[27] = (i  ,j-2,0,t6)
                outarray[28] = (i-2,j  ,0,t6)
                outarray[29] = (i-2,j+2,0,t6)
                
                #t7 hopping
                outarray[30] = (i-2,j+2,1,t7)
                outarray[31] = (i-1,j+2,1,t7)
                outarray[32] = (i+2,j-1,1,t7)
                outarray[33] = (i+2,j-2,1,t7)
                outarray[34] = (i-2,j-1,1,t7)
                outarray[35] = (i-1,j-2,1,t7)
                
                #t8 hopping
                outarray[36] = (i+1,j+1,1,t8)
                outarray[37] = (i+1,j-3,1,t8)
                outarray[38] = (i-3,j+1,1,t8)
                
            if s == 1:
                #t1 hopping
                outarray[0]  = (i  ,j  ,0,t1)
                outarray[1]  = (i+1,j  ,0,t1)
                outarray[2]  = (i  ,j+1,0,t1)
        
                #t2 hopping
                outarray[3]  = (i  ,j+1,1,t2)
                outarray[4]  = (i+1,j  ,1,t2)
                outarray[5]  = (i+1,j-1,1,t2)
                outarray[6]  = (i  ,j-1,1,t2)
                outarray[7]  = (i-1,j  ,1,t2)
                outarray[8]  = (i-1,j+1,1,t2)
                
                #t3 hopping
                outarray[9]  = (i+1,j+1,0,t3)
                outarray[10] = (i+1,j-1,0,t3)
                outarray[11] = (i-1,j+1,0,t3)
                
                #t4 hopping
                outarray[12] = (i  ,j+2,0,t4)
                outarray[13] = (i-1,j+2,0,t4)
                outarray[14] = (i+2,j  ,0,t4)
                outarray[15] = (i+2,j-1,0,t4)
                outarray[16] = (i  ,j-1,0,t4)
                outarray[17] = (i-1,j  ,0,t4)
                
                #t5 hopping
                outarray[18] = (i-1,j+2,1,t5)
                outarray[19] = (i+1,j+1,1,t5)
                outarray[20] = (i+2,j-1,1,t5)
                outarray[21] = (i+1,j-2,1,t5)
                outarray[22] = (i-1,j-1,1,t5)
                outarray[23] = (i-2,j+1,1,t5)
                
                #t6 hopping
                outarray[24] = (i  ,j+2,1,t6)
                outarray[25] = (i+2,j  ,1,t6)
                outarray[26] = (i+2,j-2,1,t6)
                outarray[27] = (i  ,j-2,1,t6)
                outarray[28] = (i-2,j  ,1,t6)
                outarray[29] = (i-2,j+2,1,t6)
                
                #t7 hopping
                outarray[30] = (i+1,j+2,0,t7)
                outarray[31] = (i+2,j+1,0,t7)
                outarray[32] = (i+2,j-2,0,t7)
                outarray[33] = (i+1,j-2,0,t7)
                outarray[34] = (i-2,j+1,0,t7)
                outarray[35] = (i-2,j+2,0,t7)
                
                #t8 hopping
                outarray[36] = (i-1,j+3,0,t8)
                outarray[37] = (i+3,j-1,0,t8)
                outarray[38] = (i-1,j-1,0,t8)
            
            return outarray
                
            
            '''
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
                                
                                #need a control statement to determine if its t3 or t4
                                out_array[k] = [i3,j3,s3,t3]
                                out_array[k] = [i3,j3,s3,t4]
                                k += 1
            return out_array     
            '''
            
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
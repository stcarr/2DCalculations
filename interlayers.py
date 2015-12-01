# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 17:31:28 2015

@author: Stephen
"""
import math

class Interlayer:
        
        # Initialising Interlayer does not change any of its properties
        # this is used as a class which returns interlayer hamiltonian energies
        # and so it is used as a sort of library look-up by heterostructure.py
        def __init__(self):
            null = 0
            
        def interaction(self,mat1,mat2,grid1,grid2,r,a1,a2,r12,r21):
            if mat1 == 'graphene':
                if mat2 == 'graphene':
                    return self.int_graphene_graphene(grid1,grid2,r,a1,a2,r12,r21)
            
        def int_graphene_graphene(self,grid1,grid2,r,a1,a2,r12,r21):
        
            a = 2.46

            lambda0 = 0.3155
            lambda3 = -0.0688
            lambda6 = -0.0083
            
            if r == 0:
                return lambda0
            
            xi0 = 1.7543
            xi3 = 3.4692
            xi6 = 2.8764
            
            x0 = False
            x3 = 0.5212
            x6 = 1.5206
            
            kappa0 = 2.0010
            kappa3 = False
            kappa6 = 1.5731
            
            pre12 = (r**2+a1**2-r12**2)/(2*r*a1)
            pre21 = (r**2+a2**2-r21**2)/(2*r*a2)
            
            if pre12 >= 1.0:
                theta12 = 0
            elif pre12 <= -1.0:
                theta12 = math.pi
            else:
                theta12 = math.acos((r**2+a1**2-r12**2)/(2*r*a1))
                
            if pre21 >= 1.0:
                theta21 = 0
            elif pre21 <= -1.0:
                theta21 = math.pi
            else:
                theta21 = math.acos((r**2+a2**2-r21**2)/(2*r*a2))

            
            def V0(r):
                return lambda0*math.exp(-xi0*(r/a)**2)*math.cos(kappa0*(r/a))
            
            def V3(r):
                return lambda3*((r/a)**2)*math.exp(-xi3*((r/a)-x3)**2)
                
            def V6(r):
                return lambda6*math.exp(-xi6*((r/a)-x6)**2)*math.sin(kappa6*(r/a))
            
            energy = V0(r) + V3(r)*(math.cos(3*theta12) + math.cos(3*theta21)) + V6(r)*(math.cos(6*theta12) + math.cos(6*theta21))
        
            return energy    
                
if __name__ == '__main__':
    ilay = Intralayer()
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 15:35:56 2015

@author: Stephen
"""

class Locality:
    
    def __init__(self):
        self.vals = []
        self.vecs = []
        
    def addVals(self,vals):
        self.vals.append(vals)
    
    def addVecs(self,vecs):
        self.vecs.append(vecs)
        
    def getVals(self):
        return self.vals
        
    def getVecs(self):
        return self.vecs
        
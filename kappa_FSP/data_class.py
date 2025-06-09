
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 11:31:57 2024

@author: alfie
"""

import numpy as np
import math

class data:
    
    def __init__(self, n_input, r_input, k_input, width_input, seed_input):
        
        self.n = n_input
        self.width = width_input
        self.seed = seed_input
        self.r = r_input
        self.k = k_input
        
        # index of vertices
        self.V = range(self.n)
        self.R = range(self.n,self.n + self.r)
        self.T = range(self.r + self.n)
        
        
        # set of tree arcs
        self.A_t = [(i,j) for i in self.T for j in self.V]
        
        # set of assignment arcs
        self.A_a = [(i,j) for i in self.V for j in self.V]
        
        # set of root arcs
        self.A_r = [(i,j) for i in self.V for j in self.R]
        
        self.loc = None
        self.t = None
        self.a = None
        
    def create_data(self):
        
        rnd = np.random.RandomState(self.seed)
        
        self.loc = {i:(rnd.random()*self.width,rnd.random()*self.width) for i in self.T}
        
        self.t = {(i,j): math.hypot(self.loc[i][0]-self.loc[j][0],self.loc[i][1]-self.loc[j][1]) for (i,j) in self.A_t}
        self.a = {(i,j): 2*math.hypot(self.loc[i][0]-self.loc[j][0],self.loc[i][1]-self.loc[j][1]) for (i,j) in self.A_a}
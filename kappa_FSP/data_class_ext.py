#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 12:42:48 2025

@author: alfie
"""

import numpy as np
import math

class data_ext:
    
    def __init__(self, n_input, r_input, k_input, width_input, seed_input):
        
        self.n = n_input
        self.width = width_input
        self.seed = seed_input
        self.r = r_input
        self.k = k_input
        
        self.V = range(self.n)
        self.R = range(self.n,self.n + self.r)
        self.T = range(self.r + self.n)
        
        self.A_t = [(i,j) for i in self.T for j in self.T]
        self.A_r_t = [(r,(i,j)) for r in self.R for (i,j) in self.A_t]
        self.A_a = [(i,j) for i in self.V for j in self.V]
        self.A_r = [(i,j) for i in self.V for j in self.R]
        
        self.loc = None
        self.t = None
        self.a = None
        
    def create_data(self):
        
        rnd = np.random.RandomState(self.seed)
        
        self.loc = {i:(rnd.random()*self.width,rnd.random()*self.width) for i in self.T}
    
        self.t = {(i,j): math.hypot(self.loc[i][0]-self.loc[j][0],self.loc[i][1]-self.loc[j][1]) for (i,j) in self.A_t}
        self.a = {(i,j): 2*math.hypot(self.loc[i][0]-self.loc[j][0],self.loc[i][1]-self.loc[j][1]) for (i,j) in self.A_a}
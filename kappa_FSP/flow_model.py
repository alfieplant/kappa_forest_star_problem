#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 12:53:11 2025

@author: alfie
"""


from docplex.mp.model import Model
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from helper import *

class flow_model():

    
    def __init__(self, data, solved_model):
        
        self.model = Model("flow")
        
        x_sol = {(i,j):1 for (i,j) in data.A_t if solved_model.x[i,j].solution_value > 0.1}
        y_sol = {(i,i):1 for i in data.V if solved_model.y[i,i].solution_value > 0.1}
        
        self.keys_x = list(x_sol.keys())
        self.keys_x_big = []
        self.keys_y = list(y_sol.keys())
        
        for (i,j) in self.keys_x:
            for v in data.R:
                if i < data.n:
                    self.keys_x_big.append((i,j,v))
                elif i == v:
                    self.keys_x_big.append((i,j,v))
                    
        self.keys_phi = []
        
        for (i,j) in self.keys_x:
            for v in data.R:
                if i < data.n:
                    self.keys_phi.append((j,i,v))
                elif i == v:
                    self.keys_phi.append((j,i,v))
                    
        for i in data.R:
            self.keys_phi.append((i, data.n + data.r, i))
            
        self.keys_w = [(i,v) for (i,j) in self.keys_y for v in data.R]
            
        self.phi = self.model.continuous_var_dict(self.keys_phi, name = 'phi')
            
        self.w = self.model.binary_var_dict(self.keys_w, name = 'w')
            
        self.x = self.model.binary_var_dict(self.keys_x_big, name = 'x')
      
        self.model.add_constraints(self.model.sum(self.w[j,v] for v in data.R) == data.k*y_sol[j,j] for j in data.V if (j,j) in self.keys_y)
        
        self.model.add_constraints(self.w[j,k] + self.model.sum(self.phi[i,j,k] for i in data.V if (i,j,k) in self.keys_phi) == self.model.sum(self.phi[j,i,k] for i in data.T if (j,i,k) in self.keys_phi) for j in data.V for k in data.R if (j,j) in self.keys_y)
        
        self.model.add_constraints(self.phi[v, data.n + data.r, v]  == self.model.sum(self.w[i,v] for i in data.V) for v in data.R if (i,i) in self.keys_y)
        
        self.model.add_constraints(self.phi[i,j,k] <= data.n*self.x[j,i,k] for (i,j,k) in self.keys_phi if (j,i,k) in self.keys_x_big)
        
        self.model.add_constraints(self.model.sum(self.x[i,j,k] for k in data.R if (i,j,k) in self.keys_x_big) == x_sol[i,j] for (i,j) in self.keys_x)
        
        self.model.add_constraints(self.model.sum(self.x[j,i,g] for j in data.T if (j,i,g) in self.keys_x_big) == self.w[i,g] for i in data.V for g in data.R if (i,i) in self.keys_y)
        
        self.model.add_constraints(self.x[i,l,g] <= self.model.sum(self.x[k,i,g] for k in data.T if (k,i,g) in self.keys_x_big if k != l) for i in data.V for l in data.V for g in data.R if (i,l,g) in self.keys_x_big)
        
        self.solution = self.model.solve(log_output = True)
        
        
    def plot_solution(self, data, solved_model):
        
        plt.figure()
        plt.axis('off')
        plt.title('Flow Model')
        fig = plt.gcf()
        fig.set_size_inches(8, 8)
        cmap = cm.viridis
    
        for i in data.V:
            plt.scatter(data.loc[i][0], data.loc[i][1], c='black')
            plt.annotate(i, (data.loc[i][0] + 2, data.loc[i][1]))
    
        for i in data.R:
            plt.scatter(data.loc[i][0], data.loc[i][1], c = cmap((i - data.n) / data.r), marker = 'x')
            plt.annotate(i, (data.loc[i][0] + 2, data.loc[i][1]))
    
        for (i,j,k) in self.keys_x_big:
            if self.x[i,j,k].solution_value > 0.1:

                mid_x = (data.loc[i][0] + data.loc[j][0]) / 2
                mid_y = (data.loc[i][1] + data.loc[j][1]) / 2
        
                plt.plot([data.loc[i][0], data.loc[j][0]],[data.loc[i][1], data.loc[j][1]], c=cmap((k - data.n) / data.r), lw=1)
        
                plt.arrow(data.loc[i][0], data.loc[i][1], mid_x - data.loc[i][0], mid_y - data.loc[i][1], head_width=1, head_length=1.5, fc=cmap((k - data.n) / data.r), ec=cmap((k - data.n) / data.r))
                
                
        for (i, j) in data.A_a:
            if (i != j and solved_model.y[i, j].solution_value > 0.9):
    
                plt.plot([data.loc[i][0], data.loc[j][0]],[data.loc[i][1], data.loc[j][1]], c='red', lw=1)
    
    
        plt.show()
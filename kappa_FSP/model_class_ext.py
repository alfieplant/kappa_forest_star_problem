#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 12:44:01 2025

@author: alfie
"""

from docplex.mp.model import Model
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from helper_ext import *
import time

class fs_model_ext():
    
    def __init__(self, data):
        
        
        self.model = Model("Forest_star_extended")
        
        self.x = self.model.binary_var_dict(data.A_r_t, name = 'x')
        
        self.y = self.model.binary_var_dict(data.A_a, name = 'y')
        
        self.w = self.model.binary_var_dict(data.A_r, name = 'w')
    
        self.model.minimize(self.model.sum(self.x[r,(i,j)]*data.t[i,j] for (r,(i,j)) in data.A_r_t) + self.model.sum(self.y[i,j]*data.a[i,j] for (i,j) in data.A_a))
        
        self.model.add_constraints(self.model.sum(self.y[i,j] for j in data.V) == 1 for i in data.V)
        
        self.model.add_constraints(self.y[i,i]  >= self.y[j,i] for i in data.V for j in data.V if i != j)
        
        self.model.add_constraints(self.model.sum(self.w[i,r] for r in data.R) == data.k * self.y[i,i] for i in data.V)
        
        self.model.add_constraints(self.model.sum(self.x[(r,e)] for e in get_cutset_in_root(data.A_a, [i], r)) == self.w[i,r] for r in data.R for i in data.V)
        
        self.model.add_constraints(self.model.sum(self.x[r,(j,i)] for (j,i) in get_cutset_in_root(data.A_a, [i], r) if j != k) >= self.x[r,(i,k)] for r in data.R for k in data.V for i in data.V if i != k)
        
        self.model.add_constraints(self.model.sum(self.x[r,(i,j)] + self.x[r,(j,i)] for r in data.R) <= 1 for (i,j) in data.A_t)
        
    def solve(self,log):
        self.model.solve(log_output = log)
        

    def plot_solution(self, data):
        plt.figure()
        plt.axis('off')
        plt.title('Naive Formulation')
        fig = plt.gcf()
        fig.set_size_inches(8, 8)
        cmap = cm.viridis

        for i in data.V:
            plt.scatter(data.loc[i][0], data.loc[i][1], c='black')
            plt.annotate(i, (data.loc[i][0] + 2, data.loc[i][1]))

        for i in data.R:
            plt.scatter(data.loc[i][0], data.loc[i][1], c = cmap((i - data.n) / data.r), marker = 'x')
            plt.annotate(i, (data.loc[i][0] + 2, data.loc[i][1]))

        for k in data.R:
            for (k,(i,j)) in data.A_r_t:
                if self.x[k,(i,j)].solution_value > 0.9:
    
                    mid_x = (data.loc[i][0] + data.loc[j][0]) / 2
                    mid_y = (data.loc[i][1] + data.loc[j][1]) / 2
            
                    plt.plot([data.loc[i][0], data.loc[j][0]],[data.loc[i][1], data.loc[j][1]], c=cmap((k - data.n) / data.r), lw=1)
                    plt.arrow(data.loc[i][0], data.loc[i][1], mid_x - data.loc[i][0], mid_y - data.loc[i][1], head_width=1, head_length=1.5, fc=cmap((k - data.n) / data.r), ec=cmap((k - data.n) / data.r))
                
                
        for (i, j) in data.A_a:
            if self.y[i, j].solution_value > 0.9:
                plt.plot([data.loc[i][0], data.loc[j][0]],[data.loc[i][1], data.loc[j][1]], c='red', lw=1)


        plt.show()
        
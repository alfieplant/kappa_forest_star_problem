#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 12:27:56 2025

@author: alfie
"""

from docplex.mp.model import Model
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from helper import *
from helper import get_cutset_in
import time

class fs_model():
    
    def __init__(self, data):
    
        
        self.model = Model("Forest_star")
        
        self.x = self.model.binary_var_dict(data.A_t, name = 'x')
        
        self.y = self.model.binary_var_dict(data.A_a, name = 'y')
    
        self.model.add_constraints(self.model.sum(self.y[i,j] for j in data.V) == 1 for i in data.V)
        
        self.model.add_constraints(self.y[i,i]  >= self.y[j,i] + self.x[i,j] + self.x[j,i] for i in data.V for j in data.V if i != j)
        
        self.model.add_constraints(self.model.sum(self.x[e] for e in get_cutset_in(data.A_t, [i])) == data.k * self.y[i,i]  for i in data.V)
                                   
        self.model.add_constraints(self.x[i,j] + self.x[j,i] <= 1 for (i,j) in data.A_t if (j,i) in data.A_t)
        
        self.model.add_constraints(self.model.sum(self.x[j,i] for (j,i) in get_cutset_in(data.A_t, [i]) if j != k) >= self.x[i,k] for k in data.V for i in data.V if i != k)
        
        self.model.minimize(self.model.sum(self.x[i,j]*data.t[i,j] for (i,j) in data.A_t) + self.model.sum(self.y[i,j]*data.a[i,j] for (i,j) in data.A_a))
        
        
    def solve(self,log):
        self.model.solve(log_output = log)
        
    def plot_solution(self, data):
        plt.figure()
        plt.axis('off')
        plt.title('Reduced Formulation')
        fig = plt.gcf()
        fig.set_size_inches(8, 8)
        cmap = cm.viridis

        for i in data.V:
            plt.scatter(data.loc[i][0], data.loc[i][1], c='black')
            plt.annotate(i, (data.loc[i][0] + 2, data.loc[i][1]))

        for i in data.R:
            plt.scatter(data.loc[i][0], data.loc[i][1], c = 'black', marker = 'x')
            plt.annotate(i, (data.loc[i][0] + 2, data.loc[i][1]))

        for (i,j) in data.A_t:
            if self.x[i,j].solution_value > 0.9:

                mid_x = (data.loc[i][0] + data.loc[j][0]) / 2
                mid_y = (data.loc[i][1] + data.loc[j][1]) / 2
        
                plt.plot([data.loc[i][0], data.loc[j][0]],[data.loc[i][1], data.loc[j][1]], c='black', lw=1)
                plt.arrow(data.loc[i][0], data.loc[i][1], mid_x - data.loc[i][0], mid_y - data.loc[i][1], head_width=1, head_length=1.5, fc='black', ec='black')
                
                
        for (i, j) in data.A_a:
            if self.y[i, j].solution_value > 0.9:

                plt.plot([data.loc[i][0], data.loc[j][0]],[data.loc[i][1], data.loc[j][1]], c='red', lw=1)


        plt.show()
        
        
    def get_sol(self, data):
        
        x = []
        v_t = []
        y= []
        
        for (i,j) in data.A_t:
            if self.x[i,j].solution_value > 10e-6:
                x.append((i,j))
                
        for i in data.V:
            if self.y[i,i].solution_value > 10e-6:
                v_t.append(i)
                
        for (i,j) in data.A_a:
            if self.y[i,j].solution_value > 10e-6:
                y.append((i,j))
                
        return x, v_t, y
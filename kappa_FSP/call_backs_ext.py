#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 12:45:48 2025

@author: alfie
"""

from cplex.callbacks import *
from docplex.mp.callbacks.cb_mixin import *
from helper_ext import *
import igraph as ig
from model_class_ext import fs_model_ext
import matplotlib.pyplot as plt
import time

class Callback_lazy_ext(ConstraintCallbackMixin,LazyConstraintCallback):
    
    def __init__(self,env):
        

        LazyConstraintCallback.__init__(self,env)
        ConstraintCallbackMixin.__init__(self)
        

                 
    def __call__(self):
        
        # LP Solver
        sol_x = self.make_solution_from_vars(self.model_instance.x.values())
        sol_y = self.make_solution_from_vars(self.model_instance.y.values())
        sol_w = self.make_solution_from_vars(self.model_instance.w.values())
    
        arcs_r_t = {}
        x_vals_r = {}
        v_r_t = {}
        
        arcs_a = []
        y_vals = []
                     
        
        for (i,j) in self.data_instance_ext.A_a:
            if sol_y.get_value(self.model_instance.y[i,j]) > 10e-6:
                arcs_a.append((i,j))
                y_vals.append(min(sol_y.get_value(self.model_instance.y[i,j]),1))
                
        for k in self.data_instance_ext.R:
            
            arcs_r_t[k] = []
            x_vals_r[k] = []
            v_r_t[k] = []
            
            
            for i in self.data_instance_ext.V:
                if sol_w.get_value(self.model_instance.w[i,k]) > 10e-6:
                    v_r_t[k].append(i)

            for (i,j) in self.data_instance_ext.A_a:
                if sol_x.get_value(self.model_instance.x[k,(i, j)]) > 10e-6:
    
                
                    arcs_r_t[k].append((i,j))
                    x_vals_r[k].append(min(sol_x.get_value(self.model_instance.x[k,(i, j)]),1))
                    
            for i in self.data_instance_ext.V:            
                if sol_x.get_value(self.model_instance.x[k,(k, i)]) > 10e-6:
                    arcs_r_t[k].append((k,i))
                    x_vals_r[k].append(min(sol_x.get_value(self.model_instance.x[k,(k, i)]),1))
                    
        
            g = graph(v_r_t[k] + [k], arcs_r_t[k])
            
            components = list(g.connected_components(mode='weak'))
            
            # Check connectivity
            if len(components) > 1:
                
                for component in components:
                    
                    if max(g.vs[component]['name']) < self.data_instance_ext.n:

                        
                        cut_set = get_cutset_in_root(self.data_instance_ext.A_a, g.vs[component]['name'], k)
            
                        for i in g.vs[component]['name']:  
                            ct = self.model_instance.model.sum(self.model_instance.x[k,e] for e in cut_set) >= self.model_instance.w[i,k]
                            ct_cpx = self.linear_ct_to_cplex(ct)
                            self.add(ct_cpx[0],ct_cpx[1],ct_cpx[2])
        
                    
        
        
class Callback_user_ext(ConstraintCallbackMixin, UserCutCallback):
    
    def __init__(self,env):
        

        UserCutCallback.__init__(self,env)
        ConstraintCallbackMixin.__init__(self)
        
    def __call__(self):
        
        # LP Solver
        sol_x = self.make_solution_from_vars(self.model_instance.x.values())
        sol_y = self.make_solution_from_vars(self.model_instance.y.values())
        sol_w = self.make_solution_from_vars(self.model_instance.w.values())
    
        arcs_r_t = {}
        x_vals_r = {}
        v_r_t = {}
        
        arcs_a = []
        y_vals = []
        
        for (i,j) in self.data_instance_ext.A_a:
            if sol_y.get_value(self.model_instance.y[i,j]) > 1e-8:
                arcs_a.append((i,j))
                y_vals.append(min(sol_y.get_value(self.model_instance.y[i,j]),1))
                

        # Separation Procedure
        for v in self.data_instance_ext.R:
        
            
            arcs_r_t[v] = []
            x_vals_r[v] = []
            v_r_t[v] = []
            

            
            for i in self.data_instance_ext.V:
                if sol_w.get_value(self.model_instance.w[i,v]) > 1e-8:
                    v_r_t[v].append(i)

            for (i,j) in self.data_instance_ext.A_a:
                if sol_x.get_value(self.model_instance.x[v,(i, j)]) > 1e-6:
                    
                    arcs_r_t[v].append((i,j))
                    x_vals_r[v].append(min(sol_x.get_value(self.model_instance.x[v,(i, j)]),1))
                    
            for i in self.data_instance_ext.V:
                if sol_x.get_value(self.model_instance.x[v,(v, i)]) > 1e-6:
                    arcs_r_t[v].append((v,i))
                    x_vals_r[v].append(min(sol_x.get_value(self.model_instance.x[v,(v, i)]),1))
                    
            
            g = graph(v_r_t[v] + [v], arcs_r_t[v], x_vals_r[v])
            
            components = list(g.connected_components(mode='weak'))
            
            subtour = False
            
            if len(components) > 1:
                for component in components:
                    
                    if max(g.vs[component]['name']) < self.data_instance_ext.n:                        
                        subtour = True
                        cut_set = get_cutset_in_root(self.data_instance_ext.A_a, g.vs[component]['name'], v)
            
                        for i in g.vs[component]['name']:  
                            ct = self.model_instance.model.sum(self.model_instance.x[v,e] for e in cut_set) >= self.model_instance.w[i,v]
                            ct_cpx = self.linear_ct_to_cplex(ct)
                            self.add(ct_cpx[0],ct_cpx[1],ct_cpx[2])

                        
            if subtour == False:
                for k in g.vs:
                    
                    if k['name'] < self.data_instance_ext.n:
                        cut = g.mincut(source = g.vs[-1].index, target = k.index, capacity='capacity')
                        true_cut = sum(sol_x.get_value(self.model_instance.x[v,(i,j)]) for (i,j) in get_cutset_in_root(self.data_instance_ext.A_a, g.vs[cut.partition[1]]['name'], v))


                        if true_cut < min(sol_w.get_value(self.model_instance.w[k["name"],v]),1) - 10e-6:
                            ct = self.model_instance.model.sum(self.model_instance.x[v,e] for e in get_cutset_in_root(self.data_instance_ext.A_a, g.vs[cut.partition[1]]['name'], v)) >= sol_w.get_value(self.model_instance.w[k["name"],v])
                            ct_cpx = self.linear_ct_to_cplex(ct)
                            self.add(ct_cpx[0],ct_cpx[1],ct_cpx[2])
  
                            
                            
                            

                        
                        

            
                            
 
                
                        

        
        
        
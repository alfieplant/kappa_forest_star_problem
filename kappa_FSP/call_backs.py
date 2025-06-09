#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 12:35:07 2025

@author: alfie
"""

from cplex.callbacks import *
from docplex.mp.callbacks.cb_mixin import *
from helper import *
import igraph as ig
from model_class import fs_model
import matplotlib.pyplot as plt
import time


class Callback_lazy(ConstraintCallbackMixin,LazyConstraintCallback):
    
    def __init__(self,env):
        

        LazyConstraintCallback.__init__(self,env)
        ConstraintCallbackMixin.__init__(self)
        

    def __call__(self):
        
            
        # LP Solver
        sol_x = self.make_solution_from_vars(self.model_instance.x.values())
        sol_y = self.make_solution_from_vars(self.model_instance.y.values())
        
        arcs_t = []
        x_vals = []
        arcs_a = []
        y_vals = []
        v_t = []
        
        for (i,j) in (self.data_instance.A_t):
            if sol_x.get_value(self.model_instance.x[i, j]) > 10e-6:
                arcs_t.append((i,j))
                x_vals.append(min(sol_x.get_value(self.model_instance.x[i, j]),1))
                
        for (i,j) in (self.data_instance.A_a):
            if sol_y.get_value(self.model_instance.y[i,j]) > 10e-6:
                arcs_a.append((i,j))
                y_vals.append(min(sol_y.get_value(self.model_instance.y[i,j]),1))
                
                if i == j and sol_y.get_value(self.model_instance.y[i,j]) > 10e-6:
                    v_t.append(i)
                    
        
        
        # Separation Procedure
        
        g = graph(v_t + list(self.data_instance.R), arcs_t, x_vals)
        
    
        dummy_arcs = []
        dummy_arcs_vals = []
    
        for k in self.data_instance.R:
            dummy_arcs.append((self.data_instance.n + self.data_instance.r, k))
            dummy_arcs_vals.append(1.0)
            
        g_ST = add_arcs(g, [self.data_instance.n + self.data_instance.r], dummy_arcs, dummy_arcs_vals)
        
        for v in g_ST.vs[:len(v_t)]:
            i = v.index
            cut = g_ST.mincut(source = g_ST.vs[-1].index, target = i, capacity='capacity')
            
            if cut.value < self.data_instance.k:
                S = g.vs[cut.partition[1]]['name']
                S_r = [k for k in S if k >= self.data_instance.n]
                
                ct = self.model_instance.model.sum(self.model_instance.x[e] for e in get_cutset_in(self.data_instance.A_t, S)) >= self.data_instance.k*self.model_instance.y[g.vs[i]['name'],g.vs[i]['name']] - len(S_r)
                ct_cpx = self.linear_ct_to_cplex(ct)
                self.add(ct_cpx[0],ct_cpx[1],ct_cpx[2])
            
  
class Callback_user(ConstraintCallbackMixin, UserCutCallback):
    
    def __init__(self,env):
        

        UserCutCallback.__init__(self,env)
        ConstraintCallbackMixin.__init__(self)
        
    def __call__(self):
  
        # LP Solver
        sol_x = self.make_solution_from_vars(self.model_instance.x.values())
        sol_y = self.make_solution_from_vars(self.model_instance.y.values())
        
        arcs_t = []
        x_vals = []
        arcs_a = []
        y_vals = []
        v_t = []
        
        
        for (i,j) in (self.data_instance.A_t):
            if sol_x.get_value(self.model_instance.x[i, j]) > 10e-6:
                arcs_t.append((i,j))
                x_vals.append(min(sol_x.get_value(self.model_instance.x[i, j]),1))
                
        for (i,j) in (self.data_instance.A_a):
            if sol_y.get_value(self.model_instance.y[i,j]) > 10e-8:
                arcs_a.append((i,j))
                y_vals.append(min(sol_y.get_value(self.model_instance.y[i,j]),1))
                
                if i == j and sol_y.get_value(self.model_instance.y[i,j]) > 10e-8:
                    v_t.append(i)
        

        # Separation Procedure 
        
        g = graph(v_t + list(self.data_instance.R), arcs_t, x_vals)


        dummy_arcs = []
        dummy_arcs_vals = []
        
        for k in self.data_instance.R:
            dummy_arcs.append((self.data_instance.n + self.data_instance.r, k))
            dummy_arcs_vals.append(1.0)
        
        g_ST = add_arcs(g, [self.data_instance.n + self.data_instance.r], dummy_arcs, dummy_arcs_vals)

        for v in g_ST.vs[:len(v_t)]:
            i = v.index
            cut = g_ST.mincut(source = g_ST.vs[-1].index, target = i, capacity='capacity')
            S = g.vs[cut.partition[1]]['name']
            S_r = [k for k in S if k >= self.data_instance.n]

            if cut.value < self.data_instance.k * min(1,sol_y.get_value(self.model_instance.y[g.vs[i]['name'],g.vs[i]['name']])) - 10e-6 and len(cut.partition[1]) > 1:

                true_cut = sum(sol_x.get_value(self.model_instance.x[i,j]) for (i,j) in get_cutset_in(self.data_instance.A_t, g_ST.vs[cut.partition[1]]['name']))
 
                true_rhs = self.data_instance.k * min(1,sol_y.get_value(self.model_instance.y[g.vs[i]['name'],g.vs[i]['name']])) - len(S_r)

                ct = self.model_instance.model.sum(self.model_instance.x[e] for e in get_cutset_in(self.data_instance.A_t, S)) >= self.data_instance.k * max(1,sol_y.get_value(self.model_instance.y[g.vs[i]['name'],g.vs[i]['name']])) - len(S_r)
                ct_cpx = self.linear_ct_to_cplex(ct)
                self.add(ct_cpx[0],ct_cpx[1],ct_cpx[2])
                
                        

            
            
            
            
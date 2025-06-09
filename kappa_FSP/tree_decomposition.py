#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 12:56:32 2025

@author: alfie
"""

import pandas as pd
from data_class import data
from helper import *
import matplotlib.pyplot as plt
from model_class import *
from call_backs import *
import random
import matplotlib.cm as cm
import igraph as ig
from itertools import combinations

def tree_decomposition(data, sol):
    x_arcs = sol[0]
    v_t = sol[1].copy()

    ##### STEP 1: ASSIGN ROOT ARCS #############################################
    
    root_colours = {}
    assigned = {i : [] for i in v_t}
    
    for i, root in enumerate(data.R):
        root_colours[root] = i
        
    arc_colours = {(i,j) : root_colours[i] for (i,j) in x_arcs if i in data.R}
    
    for (i,j) in arc_colours.keys():
        assigned[j].append(root_colours[i])
        
    ##### STEP 2: ASSIGN CHOKEPOINT PATHS ######################################
    
    chokepoints_r, chokepoints_k = find_all_chokepoints(data, sol)
    print(chokepoints_r)
    
    paths = {}
    
    for k in v_t:
        paths[k] = sorted(dfs(k, list(data.R), x_arcs, path_list=[]), key=lambda x: len(x))           
            
    ##### STEP 3a: ASSIGN CHOKEPOINT PATHS #######################################
                        
                        
    for k in v_t:
        if len(assigned[k]) < data.k:
            
            for path in paths[k]:
                
                r = path[-1][0]
                if root_colours[r] not in assigned[k]:
                
                    if (k,r) in chokepoints_k.keys():
                        if all(C in path for C in chokepoints_k[k,r]):
                            assigned_arcs = [arc for arc in path if arc in arc_colours.keys()]
                            
                            if all(arc_colours[arc] != root_colours[j] for arc in assigned_arcs for j in data.R if j != r):
                                if all(C not in path for C in chokepoints_r.keys() if chokepoints_r[C] != r):
                                
                                    for arc in path:
                                        arc_colours[arc] = root_colours[r]
                                        
                                        if root_colours[r] not in assigned[arc[1]]:
                                            assigned[arc[1]].append(root_colours[r])  
           
        ##### STEP 3b: ASSIGN REMAINING PATHS #######################################
            
        if len(assigned[k]) < data.k:
            for path in paths[k]:
                r = path[-1][0]
                if root_colours[r] not in assigned[k]:
                    assigned_arcs = [arc for arc in path if arc in arc_colours.keys()]
                    
                
                    if all(arc_colours[arc] != root_colours[j] for arc in assigned_arcs for j in data.R if j != r):
                        
                        if all(C not in path for C in chokepoints_r.keys() if chokepoints_r[C] != r):
                            
                            
                                            
                            for arc in path:
                                arc_colours[arc] = root_colours[r]
                                
                                if root_colours[r] not in assigned[arc[1]]:
                                    assigned[arc[1]].append(root_colours[r])
                                    
            
                                    
        plot_paths(data, x_arcs, arc_colours, root_colours, v_t)

                
def find_all_chokepoints(data, sol):
    x_arcs = [(i,j) for (i,j) in sol[0] if i not in data.R]
    root_arcs = [(i,j) for (i,j) in sol[0] if i in data.R]
    v_t = sol[1].copy()
    
    chokepoints_r = {}
    chokepoints_k = {}
    completed_paths = {r: [] for r in data.R}
    
    arc_sets = sorted(subsets(x_arcs, 1, data.r), key=lambda x: len(x))
    potential_chokepoints = [arc_set for arc_set in arc_sets if all(arc_set[i][1] != arc_set[j][1] for i in range(len(arc_set)) for j in range(len(arc_set)) if i !=j)]


    for C in potential_chokepoints:
        for k in v_t:
                
            root_sets = subsets(list(data.R), len(C), len(C))
    
            for D in root_sets:
                graph_arcs = [(i,j) for (i,j) in sol[0] if (i,j) not in C] + [(data.n + data.r, r) for r in D]
                
                g = graph(v_t + list(range(data.n, data.n + data.r + 1)), graph_arcs, weights = [1]*len(graph_arcs))
                
                cut = g.mincut(source = g.vs[-1], target = g.vs.find(name=k).index, capacity='capacity')
                
                if cut.value == 0:

                    assigned_chokepoints = [arc for arc in C if arc in chokepoints_r.keys()]
                    unassigned_chokepoints = [arc for arc in C if arc not in chokepoints_r.keys()]
                    
                    if len(assigned_chokepoints) == 0:
                        
                        pred_assigned = {arc: [] for arc in C}
                        
                        for arc in C:
                            predecessors = [(i,arc[0]) for (i,j) in x_arcs + root_arcs if j == arc[0]]
                            for pred in predecessors:
                                if pred in chokepoints_r.keys():
                                    if chokepoints_r[pred] in D:
                                        pred_assigned[arc].append(chokepoints_r[pred])
                                        
                                if pred in root_arcs:
                                    if pred[0] in D:
                                        pred_assigned[arc].append(pred[0])
                                        
                            
                        assigned = []
                                

                        for i in range(len(D)+1):
                            for arc in C:
                                if arc not in chokepoints_r.keys():
                                    if len(pred_assigned[arc]) == i:
                                        
                                        if i == 0:
                                            for root in D:
                                                if arc not in chokepoints_r.keys() and root not in assigned:
                                                    chokepoints_r[arc] = root
                                                    assigned.append(root)
                                                    
                                                    if (k,root) in chokepoints_k.keys():
                                                        chokepoints_k[(k,root)].append(arc)
                                                    else:
                                                        chokepoints_k[(k, root)] = [arc]
                                                    
                                        else:
                                            for root in pred_assigned[arc]:
                                                if arc not in chokepoints_r.keys():
                                                    if root not in assigned:
                                                        chokepoints_r[arc] = root
                                                        assigned.append(root)
                                                    
                                                    if (k,root) in chokepoints_k.keys():
                                                        chokepoints_k[(k,root)].append(arc)
                                                    else:
                                                        chokepoints_k[(k, root)] = [arc]
                                                
                                            
                    if len(assigned_chokepoints) > 0 and len(unassigned_chokepoints) > 0:

                        roots = D.copy()
                        for arc in assigned_chokepoints:
                            roots = [r for r in roots if r != chokepoints_r[arc]]


                        for i, arc in enumerate(unassigned_chokepoints):
                            chokepoints_r[arc] = roots[i]
                                

                    if len(assigned_chokepoints) == len(C):
                        for arc in C:
                            if (k, chokepoints_r[arc]) in chokepoints_k.keys():
                                if arc not in chokepoints_k[(k, chokepoints_r[arc])]:
                                    chokepoints_k[(k,chokepoints_r[arc])].append(arc)  
                            else:
                                chokepoints_k[(k, chokepoints_r[arc])] = [arc]
                                                                              

    return chokepoints_r, chokepoints_k
    
    
def dfs(k, R, x, path_list, current_path=[], pruned=[], previous_arc=None):
    
    pruned = pruned.copy()
    current_path = current_path.copy()

    if previous_arc != None:
        current_path.append(previous_arc)
    
    predecessors = sorted([(i,k) for (i,j) in x if j == k], key=lambda x: x[0], reverse=True)
    
    if len(predecessors) == 0:
        if current_path not in pruned:
            path_list.append(current_path)
            pruned.append(current_path)
        
    
    for arc in predecessors:
        if all(arc[0] != path_arc[1] for path_arc in current_path):
            dfs(arc[0], R, x, path_list, current_path=current_path, pruned=pruned, previous_arc=arc)
        
            
        
    return path_list









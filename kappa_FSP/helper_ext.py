#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 12:45:09 2025

@author: alfie
"""

import igraph as ig
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def get_cutset_in(A, S):
    
    cutset = [(i,j) for (i,j) in A if (j in S and i not in S)]
    
    return cutset

def get_cutset_out(A, S):
    
    cutset = [(i,j) for (i,j) in A if (i in S and j not in S)]
    
    return cutset

def get_cutset_in_root(A, S, r):
    
    cutset = [(r,i) for i in S]
    cutset  += [(j,i) for (j,i) in A if (j not in S and i in S)]
    
    return cutset


def graph(vertices, arcs, weights=None):
    
    
    vertex_map = {old_label: new_index for new_index, old_label in enumerate(vertices)}
    reindexed_arcs = [(vertex_map[i], vertex_map[j]) for i, j in arcs]
    
    g = ig.Graph(n = len(vertices), edges = reindexed_arcs, edge_attrs = {'capacity': weights}, vertex_attrs = {'name': vertices}, directed=True)
    
    return g


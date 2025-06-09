#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 12:32:18 2025

@author: alfie
"""

import matplotlib.pyplot as plt
import igraph as ig
from itertools import combinations
import matplotlib.cm as cm


def get_cutset_in_root(A, S, r):
    
    cutset = [(r,i) for i in S]
    cutset  += [(j,i) for (j,i) in A if (j not in S and i in S)]
    
    return cutset


def get_cutset_in(A, S):
    
    cutset = [(i,j) for (i,j) in A if (j in S and i not in S)]
    
    return cutset

def get_cutset_out(A, S):
    
    cutset = [(i,j) for (i,j) in A if (i in S and j not in S)]
    
    return cutset


def graph(vertices, arcs, weights):
    
    vertex_map = {old_label: new_index for new_index, old_label in enumerate(vertices)}
    reindexed_arcs = [(vertex_map[i], vertex_map[j]) for i, j in arcs]
    
    g = ig.Graph(n = len(vertices), edges = reindexed_arcs, edge_attrs = {'capacity': weights}, vertex_attrs = {'name': vertices}, directed=True)
    
    return g

def add_arcs(graph, vertices, arcs, weights):
    
    n = graph.vcount()
    vertex_map = {old_label: new_index for new_index, old_label in enumerate(vertices, start=n)}
    
    for v in graph.vs:
        vertex_map[v['name']] = v.index
        
    reindexed_arcs = [(vertex_map[i], vertex_map[j]) for i, j in arcs]
    graph.add_vertices(len(vertices), attributes={'name': vertices})
    graph.add_edges(reindexed_arcs, attributes={'capacity': weights})
    
    return graph


def subsets(lst, lower, upper):
    
    subsets = []
    
    for i in range(lower,upper+1):
        subsets.extend(map(list, combinations(lst, i)))
        
    return subsets

def plot_paths(data, x, arc_colours, root_colours, v_t):   
    
    plt.figure()
    plt.axis('off')
    plt.title('Tree Decomposition')
    fig = plt.gcf()
    fig.set_size_inches(8, 8)
    cmap = cm.viridis
    
    for r in data.R:
        plt.scatter(data.loc[r][0], data.loc[r][1], c=cmap(root_colours[r] / data.r), marker = 'x')
        plt.annotate(r, (data.loc[r][0] + 2, data.loc[r][1]))
    
    for i in v_t:
        plt.scatter(data.loc[i][0], data.loc[i][1], c='black')
        plt.annotate(i, (data.loc[i][0] + 2, data.loc[i][1]))
    
    uncoloured_arcs = [(i,j) for (i,j) in x if (i,j) not in arc_colours.keys()]
    
    for (i,j) in uncoloured_arcs:
        
        mid_x = (data.loc[i][0] + data.loc[j][0]) / 2
        mid_y = (data.loc[i][1] + data.loc[j][1]) / 2
        
        plt.plot([data.loc[i][0], data.loc[j][0]],[data.loc[i][1], data.loc[j][1]], c='grey', lw=1)


        plt.arrow(data.loc[i][0], data.loc[i][1], mid_x - data.loc[i][0], mid_y - data.loc[i][1], head_width=1, head_length=1.5, fc='grey', ec='grey')
    
    
    for (i,j) in arc_colours.keys():
        
        mid_x = (data.loc[i][0] + data.loc[j][0]) / 2
        mid_y = (data.loc[i][1] + data.loc[j][1]) / 2
        
        plt.plot([data.loc[i][0], data.loc[j][0]],[data.loc[i][1], data.loc[j][1]], c=cmap(arc_colours[i,j] / data.r), lw=1)

        plt.arrow(data.loc[i][0], data.loc[i][1], mid_x - data.loc[i][0], mid_y - data.loc[i][1], head_width=1, head_length=1.5, fc=cmap(arc_colours[i,j] / data.r), ec=cmap(arc_colours[i,j] / data.r))
        
    plt.axis('off')
    plt.show()
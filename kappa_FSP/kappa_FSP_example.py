#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 13:06:53 2025

@author: alfie
"""

from data_class import data
from data_class_ext import data_ext
from model_class import *
from model_class_ext import *
from call_backs import *
from call_backs_ext import *
import pandas as pd
from tree_decomposition import *
from flow_model import *
from helper import *

# Naive Formulation

data_instance_ext = data_ext(12, 3, 2, 100, 12348) 
data_instance_ext.create_data()
naive_mdl = fs_model_ext (data_instance_ext)

cb_lazy = naive_mdl.model.register_callback(Callback_lazy_ext)
cb_lazy.model_instance = naive_mdl
cb_lazy.data_instance_ext = data_instance_ext

cb_user = naive_mdl.model.register_callback(Callback_user_ext)
cb_user.model_instance = naive_mdl
cb_user.data_instance_ext = data_instance_ext

naive_mdl.solve(True)
naive_mdl.plot_solution(data_instance_ext)

# Reduced Formulation

# data_instance = data(9, 3, 2, 100, 12431) 
# data_instance.create_data()
# reduced_mdl = fs_model(data_instance)

# cb_lazy = reduced_mdl.model.register_callback(Callback_lazy)
# cb_lazy.model_instance = reduced_mdl
# cb_lazy.data_instance = data_instance

# cb_user = reduced_mdl.model.register_callback(Callback_user)
# cb_user.model_instance = reduced_mdl
# cb_user.data_instance = data_instance

# reduced_mdl.solve(True)
# reduced_mdl.plot_solution(data_instance)

# # Tree Decomposition
# tree_decomposition(data_instance, reduced_mdl.get_sol(data_instance))

# Flow Model
# flow_mdl = flow_model(data_instance, reduced_mdl)
# flow_mdl.plot_solution(data_instance, reduced_mdl)


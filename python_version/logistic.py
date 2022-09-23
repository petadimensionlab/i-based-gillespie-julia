# -*- coding: utf-8 -*-
"""
Stochastic logistic model
@author: Shinji Nakaoka

MIT License
The MIT License (MIT)
Copyright (c) <2013> <Shinji Nakaoka>
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "individual_based_Gillespie_algorithm.cpp"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from math import exp
from ibm import *
import numpy as np

class Individual():
    def __init__(self):
        self.ID = 0
        self.group = "n"
    def i_based_reaction( self, p_id, id, group_vals_array, parm_vals, cfg ):
        r_rate = []
        r_rate.append(parm_vals[0]) # birth
        r_rate.append(parm_vals[0]*group_vals_array[0][0]/parm_vals[1]) # death
        for i in range(cfg.max_r_num):
            cfg.srv += r_rate[i]
            cfg.cr_vals[cfg.max_r_num*id+i] = cfg.srv

max_id_num = 300
max_r_id_num = 2*max_id_num

class CFG():
    def __init__(self):
        self.total_num = 0
        self.max_r_num = 0
        self.pop_num = 0
        self.id_array = np.zeros(max_id_num,dtype=int)
        self.r_id_array = np.zeros(max_r_id_num,dtype=int)
        self.cr_vals = np.zeros(max_r_id_num,dtype=float)
        self.r_types = []
        self.rid2id = {}
        self.rid2rtype = {}
        self.id2pid = {}
        self.srv = 0.0
        self.dt = 0.0
        self.r1 = 0.0
        self.r2 = 0.0
        self.selected_r_id = 0

rnd_seed = 1000
random.seed(rnd_seed)

min_step = 0.1
plot_interval = 1.0
t0 = 0.0
Tmax = 20.0

def birth_death_process( population_array, group_vals_array, group_names_array, cfg ):

    selected_r_type = cfg.rid2rtype[cfg.selected_r_id]
    selected_id = cfg.rid2id[cfg.selected_r_id]
    selected_pid = cfg.id2pid[selected_id]
    if selected_r_type=="+":
        pop_number = get_pop_number(selected_pid,group_vals_array)
        new_id = pop_number
        new_indiv = Individual()
        new_indiv.ID = pop_number
        population_array[selected_pid][new_id] = new_indiv
        update_configuration(selected_pid,selected_r_type,group_vals_array,cfg) # put update_configuration before update_number !
        update_number(group_vals_array,selected_pid,selected_r_type,"n",group_names_array)
    elif selected_r_type=="-":
        prev_group = population_array[selected_pid][selected_id].group
        update_number(group_vals_array,selected_pid,selected_r_type,prev_group,group_names_array)
        population_array[selected_pid][selected_id] = population_array[selected_pid][cfg.total_num-1]
        population_array[selected_pid][selected_id].ID = selected_id
        update_configuration(selected_pid,selected_r_type,group_vals_array,cfg)

def initialize_ibm( population_array, group_vals_array, parm_vals, cfg ):
    p_count = 0
    for pid in range(cfg.pop_num):
        m = len(group_vals_array[pid])
        g_count = 0
        for id in range(m):
            n = group_vals_array[pid][id]
            for k in range(n):
                new_indiv = Individual()
                new_indiv.ID = p_count
                new_indiv.group = "n"
                population_array[pid][k] = new_indiv
                set_initial_configuration(pid,p_count,cfg)
                p_count += 1
                g_count += 1
    cfg.total_num = p_count

cfg = CFG()
cfg.total_num = 0
cfg.max_r_num = 2
cfg.pop_num = 1
pm = ["+","-"]
for i in range(2):
    cfg.r_types.append(pm[i])

parm_num = 2
parm_vals = [1.0,100]

parm_names = {}
parm_names["r"] = 0
parm_names["K"] = 1

group_vals_array = [[]]
group_vals_array[0].append(10)

group_names_array = [{}]
group_names_array[0]["n"] = 0

population = []
for i in range(max_id_num):
    population.append("empty")
population_array = []
for i in range(cfg.pop_num):
    population_array.append(population)

initialize_ibm(population_array,group_vals_array,parm_vals,cfg)

t = t0
next_plot = plot_interval

while t<Tmax:
    #safety_check(population_array,group_vals_array,cfg)
    while next_plot<t+plot_interval:
        if next_plot>Tmax:
            break
        print "time = %f, N = %d" % (t,group_vals_array[0][0])
        next_plot += plot_interval
    i_based_Gillespie_direct(population_array,group_vals_array,parm_vals,cfg)
    birth_death_process(population_array,group_vals_array,group_names_array,cfg)
    t += cfg.dt

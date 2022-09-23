# -*- coding: utf-8 -*-
"""
Stochastic Nicholson Blowfly model
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
import csv

class Individual():
    def __init__(self):
        self.ID = 0
        self.age = 0.0
        self.group = "j"
    def i_based_reaction( self, p_id, id, group_vals_array, parm_vals, cfg ):
        r_rate = []
        if self.group=="j":
            r_rate.append( 0.0 )
            r_rate.append( parm_vals[1] ) # d_J
        elif self.group=="a":
            r_rate.append( parm_vals[0]*exp(-1.0*group_vals_array[0][1]/parm_vals[3]) )
            r_rate.append( parm_vals[2] ) # d_A
        for i in range(cfg.max_r_num):
            cfg.srv += r_rate[i]
            cfg.cr_vals[cfg.max_r_num*id+i] = cfg.srv

max_id_num = 35000
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

rnd_seed = 1
random.seed(rnd_seed)

min_step = 0.1
plot_interval = 1.0
t0 = 0.0
Tmax = 120.0

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
        update_number(group_vals_array,selected_pid,selected_r_type,"j",group_names_array)
        update_configuration(selected_pid,selected_r_type,group_vals_array,cfg)
    elif selected_r_type=="-":
        prev_group = population_array[selected_pid][selected_id].group
        update_number(group_vals_array,selected_pid,selected_r_type,prev_group,group_names_array)
        population_array[selected_pid][selected_id] = population_array[selected_pid][cfg.total_num-1]
        population_array[selected_pid][selected_id].ID = selected_id
        update_configuration(selected_pid,selected_r_type,group_vals_array,cfg)

def maturation_process( population_array, group_vals_array, parm_vals, group_names_array, dt, cfg ):
    pop_number = get_pop_number(0,group_vals_array)
    for i in range(pop_number):
        i_age = population_array[0][i].age
        i_group = population_array[0][i].group
        population_array[0][i].age += dt
        if i_age>=parm_vals[4] and i_group=="j":
            population_array[0][i].group = "a"
            update_number(group_vals_array,0,"+","a",group_names_array)
            update_number(group_vals_array,0,"-","j",group_names_array)

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
                new_indiv.group = "a"
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
#for i in range(max_id_num):
#    cfg.id_array.append(-1)
#for i in range(max_r_id_num):
#    cfg.r_id_array.append(-1)
#    cfg.cr_vals.append(0.0)

parm_num = 5
parm_vals = []
parm_vals.append(8.5) # beta
parm_vals.append(0.0060455567) # d_J
parm_vals.append(0.27) # d_A
parm_vals.append(600.0) # c
parm_vals.append(15.6) # tau

parm_names = {}
parm_names["beta"] = 0
parm_names["dJ"] = 1
parm_names["dA"] = 2
parm_names["c"] = 3
parm_names["tau"] = 4

group_vals_array = [[]]
group_vals_array[0].append(0)
group_vals_array[0].append(5000)

group_names_array = [{}]
group_names_array[0]["j"] = 0
group_names_array[0]["a"] = 1

population = []
for i in range(max_id_num):
    population.append("empty")
population_array = []
for i in range(cfg.pop_num):
    population_array.append(population)

initialize_ibm(population_array,group_vals_array,parm_vals,cfg)

t = t0
next_plot = plot_interval


# ------- debug用 (pythonスクリプトとの比較) ---------
f = open("blowfly_py.csv", "w")
writer = csv.writer(f, lineterminator='\n')
writer.writerow(["t", "J", "A"])
# --------------------------------------------------


i=0
while t<Tmax:
    
    # ------- debug用 (pythonスクリプトとの比較) ---------
    writer.writerow([t, group_vals_array[0][0], group_vals_array[0][1]])
    # --------------------------------------------------
    
    if i % 1000 == 0:
        print "time = %f, J = %d, A = %d" % (t,group_vals_array[0][0],group_vals_array[0][1])

    while next_plot<t+plot_interval:
        if next_plot>Tmax:
            break
        print "time = %f, J = %d, A = %d" % (t,group_vals_array[0][0],group_vals_array[0][1])
        count = 0
        for idv in population_array[0]:
            if not idv=="empty":
                g = idv.group
                if g=="j":
                    #print idv.ID
                    count += 1
        #print count
        next_plot += plot_interval
    i_based_Gillespie_direct(population_array,group_vals_array,parm_vals,cfg)
    birth_death_process(population_array,group_vals_array,group_names_array,cfg)
    if cfg.dt>min_step:
        dt_tmp = -cfg.dt
        t_tmp = t
        while dt_tmp<0:
            maturation_process(population_array,group_vals_array,parm_vals,group_names_array,min_step,cfg)
            dt_tmp += min_step
            t_tmp += min_step
        #resting part
        maturation_process(population_array,group_vals_array,parm_vals,group_names_array,-dt_tmp,cfg)
    else:
        maturation_process(population_array,group_vals_array,parm_vals,group_names_array,cfg.dt,cfg)    
    t += cfg.dt
    i += 1

# ------- debug用 (pythonスクリプトとの比較) ---------
f.close()
# --------------------------------------------------

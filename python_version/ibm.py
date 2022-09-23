import random, sys
from math import log, isinf

def update_number( group_vals_array, p_id, r_type, group_name, group_names_array ):
    group_id = group_names_array[p_id][group_name]
    if r_type == "+":
        group_vals_array[p_id][group_id] += 1
    elif r_type == "-":
        group_vals_array[p_id][group_id] += -1

def get_total_number( group_vals_array ):
    total_number = 0
    m = len(group_vals_array)
    for i in range(m):
        n = len(group_vals_array[i])
        for j in range(n):
            total_number += group_vals_array[i][j]
    return total_number

def get_pop_number( pid, group_vals_array ):
    pop_number = 0
    n = len(group_vals_array[pid])
    for j in range(n):
        pop_number += group_vals_array[pid][j]
    return pop_number

def get_cmr_pop_number( pid, group_vals_array ):
    pop_number = 0
    for p in range(pid+1):
        n = len(group_vals_array[p])
        for j in range(n):
            pop_number += group_vals_array[p][j]
    return pop_number
    
def safety_check( population_array, group_vals_array, cfg ):
    total_number = get_total_number(group_vals_array)
    if total_number==0:
        print "No individual is found. Stop computation."
        sys.exit()
    for pid in range(cfg.pop_num):
        pop_number = get_pop_number(pid,group_vals_array)
        if pop_number==0:
            print "No individual is found in population %s. Stop computation." % (pid)
            sys.exit()
    
    if cfg.total_num>max_id_num:
        for i in range(max_id_num):
            cfg.id_array.append(-1)
            cfg.r_id_array.append(-1)
            cfg.cr_vals.append(0.0)
        for i in range(cfg.pop_num):
            for i in range(max_id_num):
                new_indiv = Individual()
                population_array[i].append(new_indiv)

    if math.isinf(cfg.dt)==1:
        print "'dt' is infinity. Stop computation."
        sys.exit()

def pid_shift2right( pid, group_vals_array, cfg ):
    if pid==cfg.pop_num-1:
        cfg.id2pid[cfg.total_num] = pid # old: cfg.total_num+1
    else:
        for p in range(pid,cfg.pop_num-1):
            cmr_pop_number = get_cmr_pop_number(p,group_vals_array)
            cfg.id2pid[cmr_pop_number+1] = p+1
        cfg.id2pid[cfg.total_num] = cfg.pop_num-1 # old: cfg.total_num+1

def pid_shift2left( pid, group_vals_array, cfg ):
    if pid==cfg.pop_num-1:
        cfg.id2pid[cfg.total_num] = -1
    else:
        for p in range(pid,cfg.pop_num-1):
            cmr_pop_number = get_cmr_pop_number(p,group_vals_array)
            cfg.id2pid[cmr_pop_number-1] = p+1

def update_configuration( selected_pid, selected_r_type, group_vals_array, cfg ):
    if selected_r_type=="+":
        cfg.id_array[cfg.total_num] = cfg.total_num
        pid_shift2right(selected_pid,group_vals_array,cfg)
        for i in range(cfg.max_r_num):
            cfg.r_id_array[cfg.max_r_num*cfg.total_num+i] = cfg.max_r_num*cfg.id_array[cfg.total_num]+i
            cfg.rid2id[cfg.total_num*cfg.max_r_num+i] = cfg.total_num
            cfg.rid2rtype[cfg.total_num*cfg.max_r_num+i] = cfg.r_types[i]
        cfg.total_num += 1
    elif selected_r_type=="-":
        pid_shift2left(selected_pid,group_vals_array,cfg)
        cfg.total_num += -1

def set_initial_configuration( pid, p_count, cfg ):
    cfg.id_array[p_count] = p_count
    cfg.id2pid[p_count] = pid
    for i in range(cfg.max_r_num):
        cfg.r_id_array[cfg.max_r_num*p_count+i] = cfg.max_r_num*cfg.id_array[p_count]+i
        cfg.rid2id[p_count*cfg.max_r_num+i] = p_count
        cfg.rid2rtype[p_count*cfg.max_r_num+i] = cfg.r_types[i]

def generate_random_numbers( cfg ):
    cfg.r1 = random.random()
    cfg.r2 = random.random()

def i_based_Gillespie_direct( population_array, group_vals_array, parm_vals, cfg ):
    
    generate_random_numbers(cfg)
    cfg.srv = 0.0
    cmr_count = 0

    for pid in range(cfg.pop_num):
        n = len(group_vals_array[pid])
        group_count = 0
        for j in range(n):
            group_count += group_vals_array[pid][j]
        for count in range(group_count):
            population_array[pid][count].i_based_reaction(pid,cmr_count,group_vals_array,parm_vals,cfg)
            cmr_count += 1

    if cfg.srv==0.0:
        cfg.dt = min_step
        cfg.selected_r_id = -1
    else:
        thr = cfg.r2*cfg.srv
        cfg.selected_r_id = 0
        i = 1
        val = cfg.cr_vals[0]
        if val>=thr:
            cfg.selected_r_id = cfg.r_id_array[0]
        else:
            while val<thr:
                val = cfg.cr_vals[i]
                i += 1
            cfg.selected_r_id = cfg.r_id_array[i-1]
        cfg.dt = -log(cfg.r1)/cfg.srv



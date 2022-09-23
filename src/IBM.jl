module IBM

export update_number!, get_total_number, get_pop_number
export get_cmr_pop_number, pid_shift2right!, pid_shift2left!
export update_configuration!, set_initial_configuration!, i_based_Gillespie_direct
export plot_interval, t0, Tmax, max_id_num, max_r_id_num, τ, rnd_seed, min_step
export Individual, init_individual!
export CFG

include("parameters.jl")
include("Individual.jl")
include("CFG.jl")

# ------- debug用 (pythonスクリプトとの比較) ---------
export  random
using PyCall
random = pyimport("random")  # 乱数を共通化してる
# --------------------------------------------------

function update_number!(group_vals_array, p_id, r_type, group_name, group_names_array)
    group_id = group_names_array[p_id][group_name]
    if r_type == '+'
        group_vals_array[p_id][group_id] += 1
    elseif r_type == '-'
        group_vals_array[p_id][group_id] += -1
    end
end

"""
Flatten an array and return the sum
"""
function get_total_number(group_vals_array)
    return sum([(group_vals_array...)...])
end

function get_pop_number(pid, group_vals_array)
    return sum(group_vals_array[pid])
end

function get_cmr_pop_number(pid, group_vals_array)
    pop_number = 0
    for n = 1:pid+1
        pop_number += sum(group_vals_array[n])
    end
    return pop_number
end

function safety_check(population_array, group_vals_array, cfg::CFG)
    total_number = get_total_number(group_vals_array)
    if total_number == 0
        @error ("No individual is found. Stop computation.")
    end
    for pid = 1:cfg.pop_num
        pop_number = get_pop_number(pid, group_vals_array)
        if pop_number == 0
            @error ("No individual is found in population %d. Stop computation.", pid)
        end
    end

    if cfg.total_num > max_id_num
        for i = 1:max_id_num
            push!(cfg.id_array, -1)  # あっているか自信なし....
            push!(cfg.r_id_array, -1)
            push!(cfg.cr_vals, 0.0)
        end
        for _ = 1:cfg.pop_num
            for i = 1:max_id_num
                new_indiv = Individual()
                init_individual!(new_indiv)
                push!(population_array[i], new_indiv)
            end
        end
    end

    if isinf(cfg.dt)
        @error ("'dt' is infinity. Stop computation.")
    end
end

function pid_shift2right!(pid, group_vals_array, cfg::CFG)
    if pid == cfg.pop_num
        cfg.id2pid[cfg.total_num+1] = pid
    else
        for p = pid:cfg.pop_num-1
            cmf_pop_number = get_cmr_pop_number(p, group_vals_array)
            cfg.id2pid[cmr_pop_number+1] = p + 1
        end
        cfg.id2pid[cfg.total_num+2] = cfg.pop_num
    end
end

function pid_shift2left!(pid::Int64, group_vals_array, cfg::CFG)
    if pid == cfg.pop_num
        # cfg.id2pid[cfg.total_num+1] = -1
        cfg.id2pid[cfg.total_num+1] = length(group_vals_array)  # 自信なし.....
    else
        for p = pid:cfg.pop_num-1
            cmr_pop_number = get_cmr_pop_number(p, group_vals_array)
            cfg.id2pid[cmr_pop_number-1] = p + 1
        end
    end
end

function update_configuration!(selected_pid, selected_r_type, group_vals_array, cfg::CFG)
    if selected_r_type == '+'
        # shift right
        cfg.id_array[cfg.total_num+1] = cfg.total_num
        pid_shift2right!(selected_pid, group_vals_array, cfg)
        for i = 1:cfg.max_r_num
            cfg.r_id_array[cfg.max_r_num*cfg.total_num+i] =
                cfg.max_r_num * cfg.id_array[cfg.total_num+1] + i
            cfg.rid2id[cfg.total_num*cfg.max_r_num+i] = cfg.total_num + 1
            cfg.rid2rtype[cfg.total_num*cfg.max_r_num+i] = cfg.r_types[i]
        end
        cfg.total_num += 1
    elseif selected_r_type == '-'
        pid_shift2left!(selected_pid, group_vals_array, cfg)
        cfg.total_num += -1
    end
end

function set_initial_configuration!(pid, p_count, cfg::CFG)
    cfg.id_array[p_count] = p_count
    cfg.id2pid[p_count] = pid
    for i = 1:cfg.max_r_num
        cfg.r_id_array[cfg.max_r_num*(p_count-1)+i] =
            cfg.max_r_num * (cfg.id_array[p_count] - 1) + i
        cfg.rid2id[(p_count-1)*cfg.max_r_num+i] = p_count
        cfg.rid2rtype[(p_count-1)*cfg.max_r_num+i] = cfg.r_types[i]
    end
end

function generate_random_numbers!(cfg::CFG)
    # ------- debug用 (pythonスクリプトとの比較) ---------
    cfg.r1 = random.random()
    cfg.r2 = random.random()
    # -------------------------------------------------- 

    # # ------pythonスクリプトとの比較を行わない場合-------
    # cfg.r1 = rand()
    # cfg.r2 = rand()
    # # ------------------------------------------------
end


function i_based_reaction!(individual::Individual, id, group_vals_array, cfg::CFG)
    if individual.group == 'j'
        r_rate = [
            0.0,
            d_J, # d_J
        ]
    elseif individual.group == 'a'
        r_rate = [
            β * exp(-1.0 * group_vals_array[1][2] / c),
            d_A, # d_A
        ]
    end
    for i = 1:cfg.max_r_num
        cfg.srv += r_rate[i]
        cfg.cr_vals[cfg.max_r_num*id+i] = cfg.srv
    end

end


function i_based_Gillespie_direct(population_array, group_vals_array, cfg::CFG)
    generate_random_numbers!(cfg)
    cfg.srv = 0.0
    cmr_count = 0


    for pid = 1:cfg.pop_num
        group_count = sum(group_vals_array[pid])
        for count = 1:group_count
            i_based_reaction!(
                population_array[pid][count],
                cmr_count,
                group_vals_array,
                cfg,
            )
            cmr_count += 1
        end


        if cfg.srv == 0.0
            cfg.dt = min_step
            cfg.selected_r_id = -1
        else
            thr = cfg.r2 * cfg.srv
            i = 1
            val = cfg.cr_vals[1]
            while val < thr
                val = cfg.cr_vals[i+1]
                i += 1
            end
            cfg.selected_r_id = cfg.r_id_array[i]
            cfg.dt = -log(cfg.r1) / cfg.srv
        end
    end
end

# end of module
end

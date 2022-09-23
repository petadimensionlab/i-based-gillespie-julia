module blowfly

include("IBM.jl")
using .IBM
using Printf

# # ---------pythonスクリプトとの比較を行わない場合------
# using Random
# # --------------------------------------------------

# ------- debug用 (pythonスクリプトとの比較) ---------
using CSV
using DataFrames
# --------------------------------------------------

export cfg, init_cfg!

function init_cfg!(cfg::CFG)
    # initialize
    cfg.total_num = 0
    cfg.max_r_num = 0
    cfg.pop_num = 0
    cfg.id_array = zeros(max_id_num)
    cfg.r_id_array = zeros(max_r_id_num)
    cfg.cr_vals = zeros(max_r_id_num)
    cfg.srv = 0.0
    cfg.dt = 0.0
    cfg.r1 = 0.0
    cfg.r2 = 0.0
    cfg.selected_r_id = 0
    cfg.r_types = []
    cfg.id2pid = Dict()
    cfg.rid2id = Dict()
    cfg.rid2rtype = Dict()
end



function birth_death_process(
    population_array,
    group_vals_array,
    group_names_array,
    cfg::CFG,
)
    selected_r_type = cfg.rid2rtype[cfg.selected_r_id]
    selected_id = cfg.rid2id[cfg.selected_r_id]
    selected_pid = cfg.id2pid[selected_id]
    if selected_r_type == '+'
        pop_number = get_pop_number(selected_pid, group_vals_array)
        new_id = pop_number + 1
        new_indiv = Individual()
        init_individual!(new_indiv)
        new_indiv.ID = pop_number
        population_array[selected_pid][new_id] = new_indiv
        update_number!(
            group_vals_array,
            selected_pid,
            selected_r_type,
            'j',
            group_names_array,
        )
        update_configuration!(selected_pid, selected_r_type, group_vals_array, cfg)
    elseif selected_r_type == '-'
        prev_group = population_array[selected_pid][selected_id].group
        update_number!(
            group_vals_array,
            selected_pid,
            selected_r_type,
            prev_group,
            group_names_array,
        )
        population_array[selected_pid][selected_id] =
            population_array[selected_pid][cfg.total_num]
        population_array[selected_pid][selected_id].ID = selected_id
        update_configuration!(selected_pid, selected_r_type, group_vals_array, cfg)
    end
end

function maturation_process(population_array, group_vals_array, group_names_array, dt, cfg::CFG)
    pop_number = get_pop_number(1, group_vals_array)
    for i = 1:pop_number
        i_age = population_array[1][i].age
        i_group = population_array[1][i].group
        population_array[1][i].age += dt
        if i_age >= τ && i_group == 'j'
            population_array[1][i].group = 'a'
            update_number!(group_vals_array, 1, '+', 'a', group_names_array)
            update_number!(group_vals_array, 1, '-', 'j', group_names_array)
        end
    end
end

function initialize_ibm!(population_array, group_vals_array, cfg::CFG)
    p_count = 1
    for pid = 1:cfg.pop_num
        m = length(group_vals_array[pid])
        g_count = 0
        for id = 1:m
            n = group_vals_array[pid][id]
            for k = 1:n
                new_indiv = Individual()
                init_individual!(new_indiv)
                new_indiv.ID = p_count
                new_indiv.group = 'a'
                population_array[pid][k] = new_indiv
                set_initial_configuration!(pid, p_count, cfg)
                p_count += 1
                g_count += 1
            end
        end
    end
    # 1始まりで数えているので-1してる
    cfg.total_num = p_count - 1
end



function experiment()
    # set env

    cfg = CFG()
    init_cfg!(cfg)
    cfg.total_num = 0
    cfg.max_r_num = 2
    cfg.pop_num = 1
    group_names_array = [Dict('j' => 1, 'a' => 2)]

    # ------- debug(pythonスクリプトとの比較) 用 ---------
    random.seed(rnd_seed)
    # --------------------------------------------------

    # # ------pythonスクリプトとの比較を行わない場合-------
    # Random.seed!(rnd_seed)
    # # ------------------------------------------------

    pm = ['+', '-']
    cfg.r_types = pm

    population_array = [
        Array{Union{Individual, Nothing}}(nothing, max_id_num) for i = 1:cfg.pop_num
    ]

    group_vals_array = [[0, 5000]]

    initialize_ibm!(population_array, group_vals_array, cfg)
    @info initialize_ibm!

    t = t0
    next_plot = plot_interval
    # end set env

    # ------- debug(pythonスクリプトとの比較) 用 ---------
    t_array = Array{Float64}(undef, 0)
    J_array = Array{Int64}(undef, 0)
    A_array = Array{Int64}(undef, 0)
    # --------------------------------------------------


    i=1
    while t<Tmax
        if i % 1000 == 1
            @printf("tes time=%f, J=%d, A=%d\n", t, group_vals_array[1][1], group_vals_array[1][2])
        end
        
        # ------- debug用 (pythonスクリプトとの比較) ---------
        push!(t_array, t)
        push!(J_array, group_vals_array[1][1])
        push!(A_array, group_vals_array[1][2])
        # --------------------------------------------------

        while next_plot < t + plot_interval
            if next_plot > Tmax
                @info ("next_plot>Tmax")
                break
            end
            @printf(
                "tes time=%f, J=%d, A=%d\n",
                t,
                group_vals_array[1][1],
                group_vals_array[1][2]
            )
            # 
            # count = 0
            # for idv in population_array[1]
            #     if idv !== nothing
            #         g = idv.group
            #         if g == 'j'
            #             count += 1
            #         end
            #     end
            # end
            next_plot += plot_interval
        end


        i_based_Gillespie_direct(population_array, group_vals_array, cfg)
        birth_death_process(population_array, group_vals_array, group_names_array, cfg)


        if cfg.dt > min_step
            @info ("cfg.dt>min_step")
            @show cfg.dt
            @show min_step
            dt_tmp = -cfg.dt
            t_tmp = t
            @info cfg.dt
            @info dt_tmp
            @info t_tmp
            while dt_tmp < 0
                @info maturation_process
                maturation_process(
                    population_array,
                    group_vals_array,
                    group_names_array,
                    -dt_tmp,
                    cfg,
                )
                dt_tmp += min_step
                t_tmp += min_step
            end
            maturation_process(
                population_array,
                group_vals_array,
                group_names_array,
                -dt_tmp,
                cfg,
            )
        else
            maturation_process(
                population_array,
                group_vals_array,
                group_names_array,
                cfg.dt,
                cfg,
            )
        end
        t += cfg.dt
        i += 1
    end

    # ------- debug用 (pythonスクリプトとの比較) ---------
    CSV.write("blowfly_jl.csv", DataFrame(
        t=t_array, J=J_array, A=A_array
    ))
    # --------------------------------------------------
end
@time experiment()

end

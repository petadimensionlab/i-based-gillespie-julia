
mutable struct CFG
    total_num::Int64
    max_r_num::Int64
    pop_num::Int64
    id_array::Array{Int64,1}
    r_id_array::Array{Int64,1}
    cr_vals::Array{Float64,1}
    r_types::Array{Char}
    rid2id::Dict{Int64,Int64}
    rid2rtype::Dict{Int64,Char}
    id2pid::Dict{Int64,Int64}
    srv::Float64
    dt::Float64
    r1::Float64
    r2::Float64
    selected_r_id::Int64
    # Constructor
    CFG() = new()
end

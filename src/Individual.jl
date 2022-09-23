mutable struct Individual
    ID::Int64
    age::Float64
    group::Char
    Individual() = new()
end

function init_individual!(individual::Individual)
    individual.ID = 0
    individual.age = 0.0
    individual.group = 'j'
end

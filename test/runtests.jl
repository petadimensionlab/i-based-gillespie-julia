using IBMblowfly
using Test

include("../src/IBM.jl")
using .IBM
include("../src/blowfly.jl")
using .blowfly
include("../src/IBMblowfly.jl")
using .IBMblowfly
include("test_IBM.jl")
include("test_blowfly.jl")



@testset verbose = true "testibm" begin
    test_array_total = [[1, 2], [3, 4]]
    @test get_total_number(test_array_total) == 10
    test_array_pop = [[1, 2], [3, 4], [5, 6, 7]]
    @test get_pop_number(0, test_array_pop) == 3
    @test get_pop_number(2, test_array_pop) == 18
    test_array_array_cmr_pop = [[1, 2, 3], [4, 5], [6, 7]]
    @test get_cmr_pop_number(0, test_array_array_cmr_pop) == 6
    @test get_cmr_pop_number(2, test_array_array_cmr_pop) == 28
    # test pid shift to right/left settings
    test_cfg = cfg
    test_cfg.id2pid = [1, 2, 3, 4]
    test_group_vals_array = [[0, 4999]]
    test_pid = 0

    # pid_shift test

    # test shift to right
    pid_shift2right!(test_pid, test_group_vals_array, test_cfg)
    @test test_cfg.id2pid == [0, 2, 3, 4]

    # test shift to left
    test_cfg.id2pid = [1, 2, 3, 4, 5]
    test_cfg.total_num = 0

    pid_shift2left!(test_pid, test_group_vals_array, test_cfg)

    @test test_cfg.id2pid == [1, 2, 3, 4, 5]


end

@testset verbose = true "testblowfly" begin
    test_cfg = cfg
    test_cfg.total_num = 0
    test_population_array = []
    test_group_vals_array = []
    # test initialized
    @test cfg.total_num == 0
    @test cfg.max_r_num == 0
end

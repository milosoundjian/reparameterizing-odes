include("../src/hnf_column_normal.jl")
using Test

## All the examples are from the paper http://hal.inria.fr/hal-00668882

@testset "Check HNF column normal" begin
    # Example 2.4
    A = matrix(ZZ, [[8, 2, 15, 9, 11], [6, 0, 6, 2, 3]])
    H_true = matrix(ZZ, [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0]])
    V_true = matrix(ZZ, [[-1, -2, -2, -2, -1], [-3, -14, -7, -13, -7], [1, 1, 2, 1, 0], [0, 2, 0, 3, 0], [0, 1, 0, 0, 2]])
    H, V = hnf_with_normal_transform_column(A)
    @test H == H_true
    @test V == V_true

    # Example 4.4
    A = matrix(ZZ, [[6, 0, -4, 1, 3], [0, 3, 1, -4, 3]])
    H_true = matrix(ZZ, [[3, 2, 0, 0, 0], [0, 1, 0, 0, 0]])
    V_true = matrix(ZZ, [[1, 1, 2, 1, 0], [1, 0, -1, 2, 0], [1, 1, 3, 2, 1], [1, 0, 0, 2, 1], [0, 0, 0, 0, 1]])
    H, V = hnf_with_normal_transform_column(A)
    @test H == H_true
    @test V == V_true
end
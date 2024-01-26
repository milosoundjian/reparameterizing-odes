include("../src/hnf_column.jl")
using Test

## All the examples are from the paper http://hal.inria.fr/hal-00668882

@testset "Check HNF column" begin
    # Example 2.4
    A = matrix(ZZ, [[8, 2, 15, 9, 11], [6, 0, 6, 2, 3]])
    H_true = matrix(ZZ, [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0]])
    H1 = hnf_column(A)
    H2, V = hnf_with_transform_column(A)
    @test H1 == H_true
    @test H2 == H_true && A * V == H_true

    # Example 4.4
    A = matrix(ZZ, [[6, 0, -4, 1, 3], [0, 3, 1, -4, 3]])
    H_true = matrix(ZZ, [[3, 2, 0, 0, 0], [0, 1, 0, 0, 0]])
    H1 = hnf_column(A)
    H2, V = hnf_with_transform_column(A)
    @test H1 == H_true
    @test H2 == H_true && A * V == H_true

    # Example 4.7
    A = matrix(ZZ, [[-3, 1, 1, -1, 1], [0, -1, 0, -1, -2], [1, 0, 0, 1, 1]])
    H_true = matrix(ZZ, [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0]])
    H1 = hnf_column(A)
    H2, V = hnf_with_transform_column(A)
    @test H1 == H_true
    @test H2 == H_true && A * V == H_true

    # TODO: Add test cases for tall matrices
end

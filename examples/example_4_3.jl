# Example 4.3

include("../src/hnf_column_normal.jl")

# Matrix A
A = matrix(ZZ, [[2, 3]])

# Compute the Hermite normal form of A and the normal Hermite multiplier
H, V = hnf_with_normal_transform_column(A)
H
V

V_paper = matrix(ZZ, [[-1, 3], [1, -2]])

W = inv(V)
W_paper = inv(V_paper)
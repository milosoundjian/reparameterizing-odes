# Example 4.4

include("../src/hnf_column_normal.jl")

# Matrix A
A = matrix(ZZ, [[6, 0, -4, 1, 3], [0, 3, 1, -4, 3]])

# Compute the Hermite normal form of A and the normal Hermite multiplier
H, V = hnf_with_normal_transform_column(A)
H
V
W = inv(V)

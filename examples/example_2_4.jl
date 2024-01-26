# Example 2.4
using Nemo
include("../src/hnf_column_normal.jl")

# Matrix A
A = matrix(ZZ, [[8, 2, 15, 9, 11], [6, 0, 6, 2, 3]])

# Compute the Hermite normal form of A and the normal Hermite multiplier
H, U = hnf_with_normal_transform_column(A)

# Print the results
H
U
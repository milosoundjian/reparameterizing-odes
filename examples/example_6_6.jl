include("../src/hnf_column_normal.jl")

# We have the following differential equation
# t dz_1 / dt = z_1 (-2/3 + 1/3 z_1^5 z_2),
# t dz_2 / dt = z_2 (10/3 - 2/3 z_1^5 z_2 + z_1^2 z_2 / t)

# We have to calculate the matrix K (t, z_1, z_2) order
K = matrix(ZZ, [[0, 0, 2], [5, 5, 1], [1, 1, -1]])

H, U = hnf_with_transform(K)
H
# H has 1 row of 0s, so we take the last 1 row of U

A = U[end, :]

H1, V = hnf_with_normal_transform_column(A)

V
W = inv(V)
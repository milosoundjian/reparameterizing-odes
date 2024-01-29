# We have the following differential equation
# dz_1/dt = z_1(1 + z_1 z_2), dz_2/dt = z_2(1/t - z_1 z_2)

include("../src/hnf_column_normal.jl")

# Matrix A
A = matrix(ZZ, [[1, -1]])

H, V = hnf_with_normal_transform_column(A)
W = inv(V)

# How did we get A?

# We started by calculating F following example 5.2
# F = [t + t z_1 z_2, 1 - t z_1 z_2]

# We then write the matrix K with variable order (z_1, z_2) // we choose to ignore t
K = matrix(ZZ, [[1, 1], [1, 1]])
H, U = hnf_with_transform(K)

# A is the last row of U
A = U[end, :]

H1, V = hnf_with_normal_transform_column(A)
W = inv(V)
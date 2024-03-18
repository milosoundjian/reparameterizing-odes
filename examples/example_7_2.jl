include("../src/hnf_column_normal.jl")
include("../src/ode_to_matrix.jl")

using StructuralIdentifiability: parent_ring_change

ode = @ODEmodel(
    n'(t) = r * n * (1 - n / K)
)

K = ode_to_matrix(ode)

H, U = hnf_with_transform(K)

H

U

r = count_zero_rows(H)

A = U[end-r+1:end, :]

_, V = hnf_with_normal_transform_column(A)

V

W = inv(V)

V_n = V[:, r+1:end]

W_d = W[r+1:end, :]

n = size(V, 1)

V_n

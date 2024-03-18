using Nemo
using StructuralIdentifiability

using StructuralIdentifiability: parent_ring_change

include("hnf_column_normal.jl")
include("ode_to_matrix.jl")
include("proposition_6_2.jl")

function ode_reparameterizate(ode)
    # Get the ODE matrix K
    K = ode_to_matrix(ode)

    # 
    H, U = hnf_with_transform(K)

    r = count_zero_rows(H)

    A = U[end-r+1:end, :]

    _, V = hnf_with_normal_transform_column(A)

    W = inv(V)

    V_n = V[:, r+1:end]

    W_d = W[r+1:end, :]

    final_result = prop_6_2(ode, V_n, W_d)

    return final_result, W_d
end

ode = @ODEmodel(
    n'(t) = r * n * (1 - n / k)
)

ode_to_matrix(ode)

K = ode_to_matrix(ode)

H, U = hnf_with_transform(K)

A = U[end-1:end, :]

_, V = hnf_with_normal_transform_column(A)

V

W = inv(V)

get_F(ode)

ode_reparameterizate(ode)
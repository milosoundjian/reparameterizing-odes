include("hnf_column_normal.jl")
include("ode_to_matrix.jl")

# Verhulst Model of Logistic Growth
ode = @ODEmodel(
    n'(t) = r * n * (1 - n / k),
    o(t) = n(t)
)


# Determine the Matrix K
K = ode_to_matrix(ode)

K_paper = matrix(ZZ, [[1, 1], [0, -1], [1, 1], [0, 1]])

# Determine the Matrix A
H0, U = hnf_with_transform(K)
r = count_zero_rows(H)
A = U[end-r+1:end, :]

# Normal Column Hermite Multiplier of A
H1, V = hnf_with_normal_transform_column(A)
V

W = inv(V)

W_d = W[r+1:end, :]

# function full_pipeline(ode)
#     K = ode_to_matrix(ode)
#     H, U = hnf_with_transform(K)
#     r = count_zero_rows(H)
#     A = U[end-r+1:end, :]
#     H1, V = hnf_with_normal_transform_column(A)
#     Vi = V[:, 1:r]
#     Vn = V[:, r+1:end]
#     return (H, U, A, H1, V, Vi, Vn)
# end

full_pipeline(ode)



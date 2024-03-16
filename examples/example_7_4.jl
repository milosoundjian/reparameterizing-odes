include("../src/hnf_column_normal.jl")
include("../src/ode_to_matrix.jl")

using StructuralIdentifiability: parent_ring_change

ode = @ODEmodel(
    x'(t) = a - k * x(t) + h * x(t)^2 * y,
    y'(t) = b - h * x(t)^2 * y,
    # o(t) = x(t)
)

K = ode_to_matrix(ode)

H, U = hnf_with_transform(K)

H

r = count_zero_rows(H)

U

# A is the last r rows of U
A = U[end-r+1:end, :]


# ----------------------------------------

_, V = hnf_with_normal_transform_column(A)

V

# The first r columns of V
Vi = V[:, 1:r]

# The last n-r columns of V
Vn = V[:, r+1:end]

W = inv(V)

# The first r rows of W
Wu = W[1:r, :]

# The last n-r rows of W
Wd = W[r+1:end, :]

final_stuff(ode, Vn, Wd)

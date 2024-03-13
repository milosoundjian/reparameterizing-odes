include("../src/hnf_column_normal.jl")
include("../src/ode_to_matrix.jl")

ode = @ODEmodel(
    x'(t) = a - k * x(t) + h * x(t)^2 * y,
    y'(t) = b - h * x(t)^2 * y,
    o(t) = x(t)
)

ode_to_matrix(ode)

K = matrix(ZZ, [
    [0, 1, 0, 0, 0],
    [0, 0, 1, 0, 1],
    [1, 0, 0, 0, 0],
    [0, 0, 0, 1, 0],
    [1, 1, 1, 1, 1],
    [-1, 0, 1, 0, 2],
    [0, 0, 1, -1, 0]
])

# ABKH
K = matrix(ZZ, [
    [1, 0, 0, 0, 0],
    [0, 0, 0, 1, 0],

    [0, 1, 0, 0, 0],
    [0, 0, 1, 0, 1],
    
    [1, 1, 1, 1, 1],
    [-1, 0, 1, 0, 2],
    [0, 0, 1, -1, 0]
])

H, U = hnf_with_transform(K)

H

r = count_zero_rows(H)

U

# A is the last r rows of U
A = U[end-r+1:end, :]

A_paper = matrix(ZZ, [[1, 1, 1, 1, -1, 0, 0], [0, 0, -1, -3, 1, 1, 1]])
A * K

# ----------------------------------------

H1, V = hnf_with_normal_transform_column(A)

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

# Now W matches the example

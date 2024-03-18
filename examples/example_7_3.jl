include("../src/hnf_column_normal.jl")
include("../src/ode_to_matrix.jl")
include("../src")

# Lotka-Volterra Equation
ode = @ODEmodel(
    x'(t) = k_1 * a * x - k_2 * x * y,
    y'(t) = k_2 * x * y - k_3 * y,
)

# ode = @ODEmodel(
#     n'(t) = n * (r * (1 - n/K) - k * (p/(n+d))),
#     p'(t) = s * p * (1 - (h*p)/n)
# )

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

prop_6_2(ode, Vn, Wd)

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

# ----------------------------------------

# We add t to the elements of the ring for the original polynomial
new_ring_vars = vcat(ode.poly_ring.S, :t)

new_ring, vars = PolynomialRing(QQ, new_ring_vars)

# F will be an array of length n
F = [zero(new_ring) // one(new_ring) for _ in 1:new_ring.nvars]

for i in 1:length(ode.x_vars)
    F[i] = parent_ring_change(ode.x_equations[ode.x_vars[i]] // one(ode.poly_ring), new_ring)
end

F_new = [F[i] // vars[i] for i in 1:new_ring.nvars]

ode.x_equations[ode.x_vars[1]]

F_new = get_F(ode)

# ----------------------------------------

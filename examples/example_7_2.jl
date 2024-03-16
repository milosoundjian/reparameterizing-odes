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

function get_F(ode)
    new_ring_vars = vcat(ode.poly_ring.S, :t)
    new_ring, vars = PolynomialRing(QQ, new_ring_vars)

    F = [zero(new_ring) // one(new_ring) for _ in 1:new_ring.nvars]

    for i in 1:length(ode.x_vars)
        F[i] = parent_ring_change(ode.x_equations[ode.x_vars[i]] // one(ode.poly_ring), new_ring)
    end

    F[end] = one(new_ring) // one(new_ring)

    F_new = [F[i] // vars[i] for i in 1:new_ring.nvars]

    return F_new
end

dest_ring, dest_vars = PolynomialRing(QQ, ["y$i" for i in 1:(n - r)])

y = [dest_vars[i] // one(dest_ring) for i in 1:(n - r)]

function power_vector_matrix(vector, matrix)
    if length(vector) != size(matrix, 1)
        error("The length of the vector must be the same as the number of rows in the matrix")
    end

    res = [one(vector[1]) // one(vector[1]) for _ in 1:size(matrix, 2)]

    for i in eachindex(res)
        for j in eachindex(vector)
            res[i] *= vector[j]^matrix[j, i]
        end
    end

    return res
end

y_to_W_d = power_vector_matrix(y, W_d)

F = get_F(ode)

F_of_y_to_W_d = [evaluate(F[i], y_to_W_d) for i in eachindex(F)]

F_of_y_to_W_d_times_V_n = F_of_y_to_W_d * V_n

final_result = y .* F_of_y_to_W_d_times_V_n

function final_stuff(ode, V_n, W_d)
    n = size(V_n, 1)
    n_minus_r = size(V_n, 2)

    dest_ring, dest_vars = PolynomialRing(QQ, ["y$i" for i in 1:n_minus_r])
    
    F = get_F(ode)

    y = [dest_vars[i] // one(dest_ring) for i in 1:n_minus_r]

    y_to_W_d = power_vector_matrix(y, W_d)

    F_of_y_to_W_d = [evaluate(F[i], y_to_W_d) for i in eachindex(F)]

    F_of_y_to_W_d_times_V_n = F_of_y_to_W_d * V_n

    final_result = y .* F_of_y_to_W_d_times_V_n

    return final_result
end

final_stuff(ode, V_n, W_d)

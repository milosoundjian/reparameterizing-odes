using Nemo
using StructuralIdentifiability

using StructuralIdentifiability: parent_ring_change

function get_F(ode)
    # reorder ring order
    n_x_vars = length(ode.x_vars)
    n_params = length(ode.parameters)
    old_variables = vcat(ode.poly_ring.S[n_x_vars+1:end] ,ode.poly_ring.S[1:n_x_vars])

    new_ring_vars = vcat(old_variables, :t)
    new_ring, vars = PolynomialRing(QQ, new_ring_vars)

    F = [zero(new_ring) // one(new_ring) for _ in 1:new_ring.nvars]

    for i in 1:n_x_vars
        F[i+n_params] = parent_ring_change(ode.x_equations[ode.x_vars[i]] // one(ode.poly_ring), new_ring)
    end

    F[end] = one(new_ring) // one(new_ring)

    F_new = [F[i] // vars[i] for i in 1:new_ring.nvars]

    return F_new
end

# dest_ring, dest_vars = PolynomialRing(QQ, ["y$i" for i in 1:(n - r)])

# y = [dest_vars[i] // one(dest_ring) for i in 1:(n - r)]

function power_vector_matrix(vector, matrix)
    # if length(vector) != size(matrix, 1)
    #     error("The length of the vector must be the same as the number of rows in the matrix")
    # end

    res = [one(vector[1]) // one(vector[1]) for _ in 1:size(matrix, 2)]

    for i in eachindex(res)
        for j in eachindex(vector)
            res[i] *= vector[j]^matrix[j, i]
        end
    end

    return res
end

# y_to_W_d = power_vector_matrix(y, W_d)

# F = get_F(ode)

# F_of_y_to_W_d = [evaluate(F[i], y_to_W_d) for i in eachindex(F)]

# F_of_y_to_W_d_times_V_n = F_of_y_to_W_d * V_n

# final_result = y .* F_of_y_to_W_d_times_V_n

function prop_6_2(ode, V_n, W_d)
    n = size(V_n, 1)
    n_minus_r = size(V_n, 2)

    dest_ring, dest_vars = PolynomialRing(QQ, ["y$i" for i in 1:n_minus_r])
    
    F = get_F(ode)

    y = [dest_vars[i] // one(dest_ring) for i in 1:n_minus_r]

    y_to_W_d = power_vector_matrix(y, W_d)

    # This part has issues when we have output variables in the ode
    F_of_y_to_W_d = [evaluate(F[i], y_to_W_d) for i in eachindex(F)]

    F_of_y_to_W_d_times_V_n = F_of_y_to_W_d * V_n

    reparameterized_equations = y .* F_of_y_to_W_d_times_V_n

    # Dictionary  from y to reparameterized_equations
    reparameterized_equations_dict = Dict(zip(y, reparameterized_equations))

    # z is a list of the variables in the old equation
    n_x_vars = length(ode.x_vars)
    n_params = length(ode.parameters)

    n_x_vars = length(ode.x_vars)
    n_params = length(ode.parameters)
    old_variables = vcat(ode.poly_ring.S[n_x_vars+1:end] ,ode.poly_ring.S[1:n_x_vars])

    new_ring_vars = vcat(old_variables, :t)
    new_ring, vars = PolynomialRing(QQ, new_ring_vars)

    z = [vars[i] // one(new_ring) for i in 1:new_ring.nvars]

    z_to_V_n = power_vector_matrix(z, V_n)

    y_to_z = Dict(zip(y, z_to_V_n))

    return reparameterized_equations_dict, y_to_z
end

# prop_6_2(ode, V_n, W_d)

using StructuralIdentifiability
using Nemo
using StructuralIdentifiability: ODE

function vector_to_matrix(vec)
    n = length(vec)
    m = length(vec[1])
    res = zero_matrix(ZZ, n, m)
    for i in 1:n
        for j in 1:m
            res[i, j] = vec[i][j]
        end
    end
    return res
end

ode = @ODEmodel(
    n'(t) = r * n * (1 - n / k),
    # o(t) = n(t)
)

# # We then extract the equations from the ODE
x_equations = [ode.x_equations[x] for x in ode.x_vars]
# y_equations = [ode.y_equations[y] for y in ode.y_vars]

# We will need to iterate over each equation
function parse_equation(eq, ode)

    # If the equation is not fractional, turn it into a fractional equation
    if isa(eq, QQMPolyRingElem)
        eq = eq // one(parent(eq))
    end

    # Extract the numerator and denominator
    num = collect(exponent_vectors(numerator(eq)))
    den = collect(exponent_vectors(denominator(eq)))

    # Convert the numerator and denominator to a matrix
    M_num = vector_to_matrix(num)
    M_den = vector_to_matrix(den)

    # Remove the columns for the outputs
    n_x_vars = length(ode.x_vars)
    n_outputs = length(ode.y_vars)

    # Remove the n_x_vars + 1 to n_x_vars + n_outputs columns from the numerator and denominator
    # keeping the columns after that
    M_num = hcat(M_num[:, n_x_vars + n_outputs + 1:end], M_num[:, 1:n_x_vars])
    M_den = hcat(M_den[:, n_x_vars + n_outputs + 1:end], M_den[:, 1:n_x_vars])

    return (M_num, M_den)
end

parsed_x_equations = [parse_equation(eq, ode) for eq in x_equations]
# parsed_y_equations = [parse_equation(eq, ode) for eq in y_equations]

function handle_parsed_equation(num_den, index, ode)
    n_params = length(ode.parameters)
    # n_x_vars = length(ode.x_vars)

    num, den = deepcopy(num_den)

    # We will look at the first row of the denominator that we will set to 0 and change the rest accordingly
    row_to_zero = den[1, :]

    # If the first row of the denominator is 0 then we can skip this step
    if !iszero(row_to_zero)
        # We subract row_to_zero from all the rows in the numerator
        n_num = size(num, 1)
        for i in 1:n_num
            num[i, :] = num[i, :] - row_to_zero
        end

        # We subract row_to_zero from all the rows in the denominator
        n_den = size(den, 1)
        for i in 1:n_den
            den[i, :] = den[i, :] - row_to_zero
        end
    end

    # return (num, den)

    # Handling the dx/dt term (if not a output equation)
    if index != -1
        # Padding the numerator with ones in the last column (for t)
        padding = matrix(ZZ, [[1] for i in 1:size(num, 1)])
        num = hcat(num, padding)

        # subtract 1 from the n_params+indexth column of the numerator
        num[:, n_params+index] = num[:, n_params+index] - padding

        # Padding the denominator with zeros in the last column (for t)
        den = hcat(den, zero_matrix(ZZ, size(den, 1), 1))

        # Concatenate the numerator and denominator (except the first row of the denominator which is 0 and we don't need it)
        num_den = vcat(num, den[2:end, :])
    else # Handling the output equation
        # Padding the numerator with zeros in the last column (for t)
        padding = zero_matrix(ZZ, size(num, 1), 1)
        num = hcat(num, padding)

        # Padding the denominator with zeros in the last column (for t)
        den = hcat(den, zero_matrix(ZZ, size(den, 1), 1))

        # Concatenate the numerator and denominator (except the first row of the denominator which is 0 and we don't need it) 
        num_den = vcat(num, den[2:end, :])
    end

    return num_den
end

handled_x_equations = [handle_parsed_equation(parsed_x_equations[i], i, ode) for i in eachindex(parsed_x_equations)]
# handled_y_equations = [handle_parsed_equation(parsed_y_equations[i], -1) for i in eachindex(parsed_y_equations)]

function ode_to_matrix(ode)
    x_equations = [ode.x_equations[x] for x in ode.x_vars]
    y_equations = [ode.y_equations[y] for y in ode.y_vars]
    parsed_x_equations = [parse_equation(eq, ode) for eq in x_equations]
    parsed_y_equations = [parse_equation(eq, ode) for eq in y_equations]

    handled_x_equations = [handle_parsed_equation(parsed_x_equations[i], i, ode) for i in eachindex(parsed_x_equations)]
    handled_y_equations = [handle_parsed_equation(parsed_y_equations[i], -1, ode) for i in eachindex(parsed_y_equations)]


    # Concatenate the handled equations
    res = handled_x_equations[1]
    n = length(handled_x_equations)
    for i in 2:n
        res = vcat(res, handled_x_equations[i])
    end
    m = length(handled_y_equations)
    for i in 1:m
        res = vcat(res, handled_y_equations[i])
    end

    # add a row of ones at the end
    res = vcat(res, zero_matrix(ZZ, 1, size(res, 2)))
    # set the bottom right element to 1
    res[end, end] = 1

    # return transpose(res)
    return transpose(res)
end

ode_to_matrix(ode)
# We start with an ODE
using StructuralIdentifiability
using Nemo

# ode = @ODEmodel(
#     x'(t) = a * x - b * x * y,
#     y'(t) = -c * x + d * x * y,
# )

ode_verhulst = @ODEmodel(
    n'(t) = r * n * (1 - n / k),
)

# ode = @ODEmodel(
#     x'(t) = a - k * x + h * x^2 * y + 1 / (x + y),
#     y'(t) = b - h * x^2 * y,
# )

# NOTE: The order in the matrices is always x_vars, parameters, (t)

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

# FULL PROCESS

# We start with our ODE that we want to reparameterize

# Predator-prey model
ode = @ODEmodel(
    n'(t) = n * (r * (1 - n / K) - k * (p / (n + d))),
    p'(t) = s * p * (1 - h * (p / n)),
    o(t) = n * p # The observation
)

ode = @ODEmodel(
    n'(t) = r * n * (1 - n / k),
    o(t) = n(t)
)

# We then extract the equations from the ODE
x_equations = [ode.x_equations[x] for x in ode.x_vars]
y_equations = [ode.y_equations[y] for y in ode.y_vars]

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
    if n_outputs != 0
        # Remove the n_x_vars + 1 to n_x_vars + n_outputs columns from the numerator and denominator
        # keeping the columns after that
        M_num = hcat(M_num[:, 1:n_x_vars], M_num[:, n_x_vars + n_outputs + 1:end])
        M_den = hcat(M_den[:, 1:n_x_vars], M_den[:, n_x_vars + n_outputs + 1:end])
    end

    return (M_num, M_den)
end

parsed_x_equations = [parse_equation(eq, ode) for eq in x_equations]
parsed_y_equations = [parse_equation(eq, ode) for eq in y_equations]

function handle_parsed_equation(num_den, index)
    num, den = deepcopy(num_den)

    # We will look at the first row of the denominator that we will set to 0 and change the rest accordingly
    row_to_zero = den[1, :]

    # If the first row of the denominator is 0 then we can skip this step
    if !iszero(row_to_zero)
        # if !(row_to_zero == zero_matrix(ZZ, 1, length(row_to_zero)))
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

    # Handling the dx/dt term (if not a output equation)
    if index != -1
        # Padding the numerator with ones in the last column (for t)
        padding = matrix(ZZ, [[1] for i in 1:size(num, 1)])
        num = hcat(num, padding)

        # subtract 1 from the indexth column of the numerator
        num[:, index] = num[:, index] - padding

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

handled_x_equations = [handle_parsed_equation(parsed_x_equations[i], i) for i in eachindex(parsed_x_equations)]
handled_y_equations = [handle_parsed_equation(parsed_y_equations[i], -1) for i in eachindex(parsed_y_equations)]

function ode_to_matrix(ode)
    x_equations = [ode.x_equations[x] for x in ode.x_vars]
    y_equations = [ode.y_equations[y] for y in ode.y_vars]
    parsed_x_equations = [parse_equation(eq, ode) for eq in x_equations]
    parsed_y_equations = [parse_equation(eq, ode) for eq in y_equations]

    handled_x_equations = [handle_parsed_equation(parsed_x_equations[i], i) for i in eachindex(parsed_x_equations)]
    handled_y_equations = [handle_parsed_equation(parsed_y_equations[i], -1) for i in eachindex(parsed_y_equations)]


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

    # return transpose(res)
    return transpose(res)
end

ode_to_matrix(ode)

# --------- FOR LATER ------------------

system = @ODEmodel(
    EGFR'(t) = -reaction_1_k1 * EGF_EGFR(t) + reaction_1_k2 * EGF_EGFR(t) - EGFR(t) * EGFR_turnover + EGFR_turnover * pro_EGFR(t),
    pEGFR'(t) = -reaction_4_k1 * pEGFR(t) + reaction_9_k1 * EGF_EGFR(t) - pEGFR(t) * Akt(t) * reaction_2_k1 + reaction_3_k1 * pEGFR_Akt(t) + pEGFR_Akt(t) * reaction_2_k2,
    pEGFR_Akt'(t) = pEGFR(t) * Akt(t) * reaction_2_k1 - reaction_3_k1 * pEGFR_Akt(t) - pEGFR_Akt(t) * reaction_2_k2,
    Akt'(t) = pAkt(t) * reaction_7_k1 - pEGFR(t) * Akt(t) * reaction_2_k1 + pEGFR_Akt(t) * reaction_2_k2,
    pAkt'(t) = -pAkt(t) * reaction_7_k1 - pAkt(t) * reaction_5_k1 * S6(t) + reaction_6_k1 * pAkt_S6(t) + reaction_3_k1 * pEGFR_Akt(t) + pAkt_S6(t) * reaction_5_k2,
    S6'(t) = pS6(t) * reaction_8_k1 - pAkt(t) * reaction_5_k1 * S6(t) + pAkt_S6(t) * reaction_5_k2,
    pAkt_S6'(t) = pAkt(t) * reaction_5_k1 * S6(t) - reaction_6_k1 * pAkt_S6(t) - pAkt_S6(t) * reaction_5_k2,
    pS6'(t) = -pS6(t) * reaction_8_k1 + reaction_6_k1 * pAkt_S6(t),
    EGF_EGFR'(t) = reaction_1_k1 * EGF_EGFR(t) - reaction_9_k1 * EGF_EGFR(t) - reaction_1_k2 * EGF_EGFR(t),
    y1(t) = pEGFR(t) * a1 + a1 * pEGFR_Akt(t),
    y2(t) = a2 * pAkt(t) + a2 * pAkt_S6(t),
    y3(t) = pS6(t) * a3
)

ode_to_matrix(system)

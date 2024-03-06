# We start with an ODE
using StructuralIdentifiability
using Nemo

# ode = @ODEmodel(
#     x'(t) = a * x - b * x * y,
#     y'(t) = -c * x + d * x * y,
# )

# ode = @ODEmodel(
#     x'(t) = a - k * x + h * x^2 * y + 1 / (x + y),
#     y'(t) = b - h * x^2 * y,
# )

# FULL PROCESS

# We start with our ODE that we want to reparameterize

# Predator-prey model
ode = @ODEmodel(
    n'(t) = n * (r * (1 - n / K) - k * (p / (n + d))),
    p'(t) = s * p * (1 - h * (p / n))
)

# We then extract the equations from the ODE
equations = [ode.x_equations[x] for x in ode.x_vars]

# We will need to iterate over each equation
function parse_equations(equations)
    res = []
    for eq in equations
        # We determine if the equation is fractional or not
        if string(typeof(equations[1])) == "AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}"
            # If it is, we extract the numerator and denominator
            num = collect(exponent_vectors(numerator(eq)))
            den = collect(exponent_vectors(denominator(eq)))
            push!(res, (num, den))
        elseif string(typeof(equations[1])) == "QQMPolyRingElem"
            # ? You could also divide by one(ring) 
            # Else we just extract the equation
            push!(res, collect(exponent_vectors(eq)))
        else
            error("The equation is not of the correct type")
        end
    end
    return res
end

parsed_equations = parse_equations(equations)

function handle_fractional_equation(num_den)
    num, den = deepcopy(num_den)
    
    # We look at the first row of the denominator that we will set to 0 and change the rest accordingly
    row_to_zero = den[1]

    # We subract row_to_zero from all the rows in the numerator
    n_num = length(num)
    for i in 1:n_num
        num[i] = num[i] - row_to_zero
    end

    # We subract row_to_zero from all the rows in the denominator
    n_den = length(den)
    for i in 1:n_den
        den[i] = den[i] - row_to_zero
    end
    
    # We return the new numerator and denominator
    return num, den
end

num_den = handle_fractional_equation(parsed_equations[1])



# --------------------------------------------------

system = @ODEmodel(
	EGFR'(t) = -reaction_1_k1*EGF_EGFR(t) + reaction_1_k2*EGF_EGFR(t) - EGFR(t)*EGFR_turnover + EGFR_turnover*pro_EGFR(t),
	pEGFR'(t) = -reaction_4_k1*pEGFR(t) + reaction_9_k1*EGF_EGFR(t) - pEGFR(t)*Akt(t)*reaction_2_k1 + reaction_3_k1*pEGFR_Akt(t) + pEGFR_Akt(t)*reaction_2_k2,
	pEGFR_Akt'(t) = pEGFR(t)*Akt(t)*reaction_2_k1 - reaction_3_k1*pEGFR_Akt(t) - pEGFR_Akt(t)*reaction_2_k2,
	Akt'(t) = pAkt(t)*reaction_7_k1 - pEGFR(t)*Akt(t)*reaction_2_k1 + pEGFR_Akt(t)*reaction_2_k2,
	pAkt'(t) = -pAkt(t)*reaction_7_k1 - pAkt(t)*reaction_5_k1*S6(t) + reaction_6_k1*pAkt_S6(t) + reaction_3_k1*pEGFR_Akt(t) + pAkt_S6(t)*reaction_5_k2,
	S6'(t) = pS6(t)*reaction_8_k1 - pAkt(t)*reaction_5_k1*S6(t) + pAkt_S6(t)*reaction_5_k2,
	pAkt_S6'(t) = pAkt(t)*reaction_5_k1*S6(t) - reaction_6_k1*pAkt_S6(t) - pAkt_S6(t)*reaction_5_k2,
	pS6'(t) = -pS6(t)*reaction_8_k1 + reaction_6_k1*pAkt_S6(t),
	EGF_EGFR'(t) = reaction_1_k1*EGF_EGFR(t) - reaction_9_k1*EGF_EGFR(t) - reaction_1_k2*EGF_EGFR(t),
	y1(t) = pEGFR(t)*a1 + a1*pEGFR_Akt(t),
	y2(t) = a2*pAkt(t) + a2*pAkt_S6(t),
	y3(t) = pS6(t)*a3
)

equations = [system.x_equations[x] for x in system.x_vars]

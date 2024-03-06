include("hnf_column_normal.jl")
using StructuralIdentifiability
# using ModelingToolkit

# StructuralIdentifiability
ode = @ODEmodel(
    x'(t) = a * x - b * x * y,
    y'(t) = -c * x + d * x * y,
)

ode.parameters
ode.x_vars

# # ModelingToolkit
# @parameters a, b, c, d
# @variables t, x(t), y(t)

# D = Differential(t)

# eqs = [
#     D(x) ~ a * x - b * x * y,
#     D(y) ~ -c * x + d * x * y,
# ]

# measured_quantities = [o ~ x]

# ode_mtk = ODESystem(eqs, t, name = :model_name)

# ode_mtk

function ode_to_matrix()

end

function reparameterize_ode()

end

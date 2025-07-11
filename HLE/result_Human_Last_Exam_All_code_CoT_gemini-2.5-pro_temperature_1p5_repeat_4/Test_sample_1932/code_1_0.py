import sympy

# Define the symbols for the parameters
pi, d, h, rho, t, g, H = sympy.symbols('pi d h rho t g H')

# Define the expression for the mass flow rate (dm/dt)
mass_flow_rate_str = '(pi * d**2 * h * rho) / (4 * t)'

# Define the expression for the weight change (Delta W)
# This is based on the negative of the weight of the falling sand column
delta_w_expr = f"-({mass_flow_rate_str}) * sqrt(2 * g * H)"

# To display the final expression clearly, let's break it down and print it.
# dm/dt = (total mass) / (total time)
total_mass = (pi * d**2 * h * rho) / 4

# The change in weight is approximately the negative of the weight of the sand in flight.
# ΔW ≈ - (dm/dt) * v_impact_g = - (dm/dt) * g * t_fall = - (dm/dt) * g * sqrt(2H/g) = - (dm/dt) * sqrt(2gH)
# Substituting dm/dt:
# ΔW = - (π * d**2 * h * ρ / (4 * t)) * sqrt(2 * g * H)

print("The estimated change in weight, ΔW, is given by the expression:")
print("ΔW = - (π * d**2 * h * ρ / (4 * t)) * sqrt(2 * g * H)")
print("\nBreaking it down into the components of the final answer choice:")
part1 = "pi * d**2 * h * rho"
part2 = "4 * t"
part3 = "sqrt(2 * g * H)"

print(f"ΔW = - ( {part1} / ({part2}) ) * {part3}")

# To match the format of the options, we can represent it as a fraction.
print("\nOr, written as a single fraction:")
final_expression_numerator = f"pi * d**2 * h * rho * sqrt(2 * g * H)"
final_expression_denominator = "4 * t"
print(f"ΔW = - ( {final_expression_numerator} ) / ( {final_expression_denominator} )")

# Final verification against option B: -(π d² h ρ / 4t) * sqrt(2gH)
# which can be written as - (π d² h ρ * sqrt(2gH)) / (4t)
# Our result matches.
